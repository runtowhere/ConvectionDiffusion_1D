//
//  newSolver.cpp
//  ConvectionDiffusion_1D
//
//  Created by ikon on 11/6/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#include "newSolver.h"
void FirstOrderCDSolver::UpdateBoundary(mVector& u, double atTime){
    
}
void FirstOrderCDSolver::SetInitialValue(mVector& u){
    int i;
    for (i = 0; i != u.dim(); i++) {
        u[i] = mpCondition -> mpInitialFunc(0, mpCondition -> xMin + i * xStep);
    }
}

void FirstOrderCDSolver::CentralExplicitSolve(){
    
}
void FirstOrderCDSolver::ComputeCentral(mVector& uPre, mVector& uPost, double atTime){
    
}
    
void FirstOrderCDSolver::UpwindSolve(){
    SetXStep();
    mVector uPre(nNodes);
    mVector uPost(nNodes);
    SetInitialValue(uPost);
    double tNow = 0;
    while (tNow < finalTime) {
        double dt = 0.001;
        uPre = uPost;
//        ComputeTimeStep(tNow, uPre, dt);
        ComputeUpWind(uPre, uPost, tNow, dt);
        tNow += dt;
    }
    std::cout << uPost << endl;
}
void FirstOrderCDSolver::ComputeUpWind(mVector& uPre, mVector& uPost, double atTime, double& dt){
    pFunc a = mpPDE -> mFuncA;
    pFunc c = mpPDE -> mFuncC;
    double xMin = mpCondition -> xMin;
    for (int i = 1; i != nNodes - 1; i++) {
        double v = a(atTime, xMin + i * xStep, uPre[i]) * dt / xStep;
        double u = (c(atTime, xMin + i * xStep, uPre[i])
                 + std::abs(a(atTime, xMin + i * xStep, uPre[i])) * xStep / 2.0)
                 / (xStep * xStep);
        uPost[i] = (u - v / 2.0) * uPre[i + 1] + (1 - 2 * u) * uPre[i] + (u + v / 2.0) * uPre[i - 1];
    }
}

void FirstOrderCDSolver::ComputeTimeStep(double tNow, mVector& uPre, double& dt){
    double xMin = mpCondition -> xMin;
    double tMin = 0.2 * (xStep * xStep) / ( 2 * mpPDE -> mFuncC(tNow, xMin, uPre[0])
                                  + std::abs(mpPDE -> mFuncA(tNow, xMin, uPre[0])) * xStep);
    for (int i = 1; i != nNodes - 2; i++) {
        double t;
        t = 0.2 * (xStep * xStep) / ( 2 * mpPDE -> mFuncC(tNow, xMin + i * xStep, uPre[i])
                               + std::abs(mpPDE -> mFuncA(tNow, xMin + i * xStep, uPre[i])) * xStep);
        if (tMin > t) {
            tMin = t;
        }
        if (tMin < 10e-10) {
            cout << (xStep * xStep) << " h " << 2 * (mpPDE -> mFuncC(tNow, xMin + i * xStep, uPre[i]))
            << " h " <<  uPre[i] << endl;
        }
    }
    dt = tMin;
    cout << endl;
    if (finalTime - tNow < dt) {
        dt = finalTime - tNow;
    }
}

