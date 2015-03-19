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

#include <default_gui_model.h>

class ICell : public DefaultGUIModel
{
    public:
    
        ICell(void);
        virtual ~ICell(void);
        
        virtual void execute(void);
        
    protected:
    
        virtual void update(DefaultGUIModel::update_flags_t);
    
    private:
    
        //functions
        void derivs(double *,double *);
        void solve(double, double *,double *);
        
        // Model Parameters 
        static const double VNa=50.0,VK=-100.0,GNa=100.0,GK=80.0,Cm=1.0;
        static const double VL=-67.0,GL=0.1;
        static const double J_L,J_Na,J_K;
        
        double period;
        double y[4];
        double dydt[4];
        
        double V0;
        
        double Iapp_offset;
        double rate;
        int steps;
};
