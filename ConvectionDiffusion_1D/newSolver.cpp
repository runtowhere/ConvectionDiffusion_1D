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
    double tNow = initialTime;
    mVector uPre(nNodes);
    mVector uPost(nNodes);
    SetInitialValue(uPost);
    while (tNow < finalTime) {
        uPre = uPost;
        double dt = ComputeTimeStep(tNow, uPre);
        cout << dt << " t " << endl;
        ComputeUpWind(uPre, uPost, tNow, dt);
        tNow += dt;
    }
    cout << uPost << endl;
}
void FirstOrderCDSolver::ComputeUpWind(mVector& uPre, mVector& uPost, double atTime, double& dt){
    pFunc fa = mpPDE -> mFuncA;
    pFunc fc = mpPDE -> mFuncC;
    for (int i = 1; i != nNodes - 1; i++) {
        double u = (fc(atTime, xPosition(i), uPre[i]) + 0.5 * std::abs(fa(atTime, xPosition(i), uPre[i]) ) * xStep) * dt / (xStep * xStep);
        double v = fa(atTime, xPosition(i), uPre[i]) * dt / xStep;
        uPost[i] = (u - 0.5 * v) * uPre[i + 1] + (1 - 2 * u) * uPre[i] + (u + 0.5 * v) * uPre[i - 1];
    }
}

double FirstOrderCDSolver::ComputeTimeStep(double tNow, mVector& uPre){
    double controller = 0.5;
    double tMin = controller * ComputeTimeConstraint(tNow, uPre, 0);
    for (int i = 1; i != nNodes; i++) {
        double t_i = controller * ComputeTimeConstraint(tNow, uPre, i);
        if (t_i < tMin) {
            tMin = t_i;
        }
    }
    if (tNow + tMin > finalTime) {
        tMin = finalTime - tNow;
    }
    return tMin;
}

double FirstOrderCDSolver::xPosition(int i){
    return (mpCondition -> xMin + i * xStep);
}

double FirstOrderCDSolver::ComputeTimeConstraint(double tNow, mVector& uPre, int i){
    pFunc a = mpPDE -> mFuncA;
    pFunc c = mpPDE -> mFuncC;
    return (xStep * xStep) / (2 * c(tNow, xPosition(i), uPre[i]) + std::abs(a(tNow, xPosition(i), uPre[i]) ) * xStep );
}