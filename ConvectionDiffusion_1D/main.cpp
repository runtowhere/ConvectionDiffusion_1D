//
//  main.cpp
//  ConvectionDiffusion_1D
//
//  Created by ikon on 10/31/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <cassert>
#include "Matrix.h"
#include "ConvectionDiffusion_1D.h"
#include "Tridiagonal.h"

#define C1 0.1
#define Delta 0.0001

using std::cout;
using std::endl;
using std::cin;

double a(double t, double x, double u){
    return u;
}

double c(double t, double x, double u){
    return C1;
}

double rhs(double t, double x, double u){
    return 0;
}

double initial(double t, double x){
    if (x >= 0 ) {
        return 1 + Delta;
    }
    else{
        return 1;
    }
}

int main(int argc, const char * argv[])
{
    
    FirstOrderCD cd1(&a, &c, &rhs); // a, c, rhs
    BoundaryConditions_1D cond1(&initial,&initial); // boundary, initial
    cond1.SetDirichlet();
    cond1.SetRegion(-1, 1);
    FirstOrderCDSolver sol1(&cd1, &cond1);
    
    sol1.SetInitialTime(0);
    sol1.SetFinalTime(5);
    sol1.SetNumberNodes(80);
    sol1.SetTimeStep(0.001);
//    sol1.CentralExplicitSolve();
    sol1.UpwindSolve();

    
    
    return 0;
}

