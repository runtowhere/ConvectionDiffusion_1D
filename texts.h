//
//  texts.h
//  ConvectionDiffusion_1D
//
//  Created by ikon on 11/6/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#ifndef ConvectionDiffusion_1D_texts_h
#define ConvectionDiffusion_1D_texts_h

//
//  ConvectionDiffusion_1D.cpp
//  ConvectionDiffusion_1D
//  Created by Li Xinrui on 11/1/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//


void FirstOrderCDSolver::SetInitialValue(mVector &u){
    int i;
    for (i = 0; i != u.dim(); i++) {
        u[i] = mpCondition -> mpInitialFunc(0, mpCondition -> xMin + i * xStep);
    }
}

void FirstOrderCDSolver::CentralExplicitSolve(){
    SetXStep();
    mVector uPre(nNodes);
    mVector uPost(nNodes);
    SetInitialValue(uPost);
    //    cout << uPost;
    double tNow;
    for (int i = 1; (tNow = i * timeStep + initialTime) <= finalTime; i++) {
        if (ShouldChangeTimeStep(tNow, i, uPre)) {
            ;
        }
        ComputeCentral(uPre, uPost, tNow);
    }
    cout << uPost << endl;
}

void FirstOrderCDSolver::ComputeCentral(mVector& uPre, mVector& uPost, double atTime){
    double v,u;
    int i;
    uPre = uPost;
    for (i = 1; i != nNodes - 1; i++) {
        v = mpPDE -> mFuncA(atTime, mpCondition -> xMin + i * xStep, uPre[i]) * timeStep / xStep;
        u = mpPDE -> mFuncC(atTime, mpCondition -> xMin + i * xStep, uPre[i]) * timeStep / xStep * xStep;
        uPost[i] = (u - 1/2.0 * v) * uPre[i + 1] + (1 - 2 * u) * uPre[i]  + (u + 1 / 2.0 * v) * uPre[i - 1];
    }
    uPost[0] = mpCondition -> mpBoundaryFunc(atTime, mpCondition -> xMin);
    uPost[nNodes - 1] = mpCondition -> mpBoundaryFunc(atTime, mpCondition -> xMax);
}


void FirstOrderCDSolver::UpwindSolve(){
    SetXStep();
    mVector uPre(nNodes);
    mVector uPost(nNodes);
    SetInitialValue(uPost);
    double tNow;
    for (int i = 1; (tNow = i * timeStep + initialTime) <= finalTime; i++) {
        ComputeUpWind(uPre, uPost, tNow);
    }
    cout << uPost;
}
void FirstOrderCDSolver::ComputeUpWind(mVector& uPre, mVector& uPost, double atTime){
    double v,u;
    int i;
    uPre = uPost;
    for (i = 1; i != nNodes - 1; i++) {
        v = mpPDE -> mFuncA(atTime, mpCondition -> xMin + i * xStep, uPre[i]) * timeStep / xStep;
        u = (mpPDE -> mFuncC(atTime, mpCondition -> xMin + i * xStep, uPre[i])
             + std::abs(mpPDE -> mFuncA(atTime, mpCondition -> xMin + i * xStep,uPre[i]) * xStep / 2.0))
        * timeStep / xStep * xStep;
        uPost[i] = (u - 1/2.0 * v) * uPre[i + 1] + (1 - 2 * u) * uPre[i]  + (u + 1 / 2.0 * v) * uPre[i - 1];
    }
    uPost[0] = mpCondition -> mpBoundaryFunc(atTime, mpCondition -> xMin);
    uPost[nNodes - 1] = mpCondition -> mpBoundaryFunc(atTime, mpCondition -> xMax);
}

#endif
