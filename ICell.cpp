/*
 * Copyright (C) 2006 Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <math.h>
#include <ICell.h>

extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new ICell();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "Iapp",
        "A",
        DefaultGUIModel::INPUT,
    },
    {
        "V",
        "V",
        DefaultGUIModel::OUTPUT,
    },
    {
        "Iapp_offset",
        "uA/cm^2 - Current added to the input.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "rate",
        "Hz - The rate of integration.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER,
    },
    {
        "m",
        "Sodium Activation",
        DefaultGUIModel::STATE,
    },
    {
        "h",
        "Sodium Inactivation",
        DefaultGUIModel::STATE,
    },
    {
        "n",
        "Potassium Activation",
        DefaultGUIModel::STATE,
    },
};

static size_t num_vars = sizeof(vars)/sizeof(DefaultGUIModel::variable_t);

ICell::ICell(void)
    : DefaultGUIModel("ICell",::vars,::num_vars), period(RT::System::getInstance()->getPeriod()*1e-6), steps(static_cast<int>(ceil(period/25e-3))), V0(-65.0), Iapp_offset(0.0), rate(40000) {

    //V0 = -65;

    update(INIT);
    refresh();
}

ICell::~ICell(void) {}

/*******************
 * Model Functions *
 *******************/

static inline double am(double v)
{
  if(fabs(v+54)<0.001) /* use linear approx. when v is close to zero-pole */
    return 1.28+0.16*(v+54);
  else
    return 0.32*(54+v)/(1-exp(-(v+54)/4));
}

static inline double bm(double v)
{
  if(fabs(v+27)<0.001)
    return 1.4-0.14*(v+27);
  else
    return 0.28*(v+27)/(exp((v+27)/5)-1);
}

static inline double ah(double v)
{
  return 0.128*exp(-(50+v)/18);
}

static inline double bh(double v) 
{
  return 4/(1+exp(-(v+27)/5));
}

static inline double an(double v)
{
  if(fabs(v+52)<0.001)
    return 0.16+0.016*(v+52);
  else
    return .032*(v+52)/(1-exp(-(v+52)/5));
}

static inline double bn(double v)
{
  return 0.5*exp(-(57+v)/40);
}

/************************
 * Simple RK solver. *
 ************************/

 void ICell::solve(double dt, double *y, double *dydt) {
     
     // Previous Euler Solver
     // double dydt[4];
     // 
     // derivs(y,dydt);
     // 
     // for(size_t i = 0;i < 4;++i)
     // y[i] += dt*dydt[i];
     
     
     
     // number of variables
     int n = 4;
     
     //intermediate vector
     double yt[n];
     double k2[n], k3[n];
     double xh, hh, h6, h;
     int i;
     
     h = dt;        //timestep
     hh = h * 0.5;  //midpoint
     xh = h + hh;   //midpoint
     h6 = h / 6.0;  
     
     
     for (i = 0; i < n; i++) {
         yt[i] = y[i] + hh * dydt[i];
     }
     
     derivs(yt, k2);
     for (i = 0; i < n; i++) {
         yt[i] = y[i] + hh * k2[i];
     }
     
     derivs(yt, k3);
     for (i = 0; i < n; i++) {
         yt[i] = y[i] + h * k3[i];
         k3[i] += k2[i];
     }
     
     derivs(yt, k2);
     for (i = 0; i < n; i++) {
         y[i] = y[i] + h6 * (dydt[i] + k2[i] + 2.0 * k3[i]);
     }
 }

/**********************************************************
 * Macros for making the code below a little bit cleaner. *
 **********************************************************/

#define V (y[0])
#define m (y[1])
#define h (y[2])
#define n (y[3])
#define dV (dydt[0])
#define dm (dydt[1])
#define dh (dydt[2])
#define dn (dydt[3])
//#define Iapp (input(0)*1e-6 + Iapp_offset)
#define Iapp (input(0)*1e12 + Iapp_offset)

void ICell::derivs(double *y,double *dydt) {
    
    double J_L=GL*(V-VL);
    double J_Na=(GNa*pow(m,3)*h)*(V-VNa);
    double J_K=GK*pow(n,4)*(V-VK);
    
    dm=am(V)*(1-m)-bm(V)*m;
    dh=ah(V)*(1-h)-bh(V)*h;
    dn=an(V)*(1-n)-bn(V)*n;
    dV=-1.0/Cm*(J_Na+J_K+J_L-Iapp);
}

void ICell::execute(void) {

    /********************************************************************
     * Because the real-time thread may run much slower than we want to *
     *   integrate we need to run multiple interations of the solver.   *
     ********************************************************************/

     for(int i = 0;i < steps;++i) {
        derivs(y, dydt);    // added have to compute dydt first before solve
        solve(period/steps,y,dydt);
     }

    output(0) = V*1e-3;
}

void ICell::update(DefaultGUIModel::update_flags_t flag) {
    switch(flag) {
      case INIT:

          setState("m",m);
          setState("h",h);
          setState("n",n);

          setParameter("V0",V0);
          setParameter("Iapp_offset",Iapp_offset);
          setParameter("rate",rate);

          V = V0;
          // set to 0 since there is no 'inf' function given
          m = 0; //m_inf(V0);
          h = 0; //h_inf(V0);
          n = 0; //n_inf(V0);

          break;
      case MODIFY:

          V0 = getParameter("V0").toDouble();
          Iapp_offset = getParameter("Iapp_offset").toDouble();
          rate = getParameter("rate").toDouble();
          steps = static_cast<int>(ceil(period*rate/1000.0));

          V = V0;

	  // set to 0 since no 'inf' function given
          m = 0; //m_inf(V0);
          h = 0; //h_inf(V0);
          n = 0; //n_inf(V0);

          break;
      case PERIOD:
          period = RT::System::getInstance()->getPeriod()*1e-6;
          steps = static_cast<int>(ceil(period*getParameter("rate").toUInt()/1000.0));
      default:
          break;
    }
}
