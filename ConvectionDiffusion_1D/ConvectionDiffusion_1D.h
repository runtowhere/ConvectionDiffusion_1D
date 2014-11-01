//
//  ConvectionDiffusion_1D.h
//  ConvectionDiffusion_1D
//
//  Created by Li Xinrui on 11/1/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#ifndef __ConvectionDiffusion_1D__ConvectionDiffusion_1D__
#define __ConvectionDiffusion_1D__ConvectionDiffusion_1D__

#include <stdio.h>
#include "Matrix.h"

using std::cout;
using std::cin;
using std::endl;

typedef double (*pFunc)(double t, double x);

class FirstOrderCD{
    friend class FirstOrderCDSolver;
private:
    pFunc mFuncA;
    pFunc mFuncC;
    pFunc mRHSFunc;
public:
    FirstOrderCD(pFunc a, pFunc c, pFunc f){
        mFuncA = a;
        mFuncC = c;
        mRHSFunc = f;
    }
};

class BoundaryConditions_1D{
    friend class FirstOrderCDSolver;
private:
    bool isDirichelt = 0;
    pFunc mpBoundaryFunc;
    pFunc mpInitialFunc;
    double xMin;
    double xMax;
public:
    void SetDirichlet(){
        isDirichelt = 1;
    }
    void SetBoundaryFunc(pFunc f){
        mpBoundaryFunc = f;
    }
    void SetInitialFunc(pFunc f){
        mpInitialFunc = f;
    }
    void SetRegion(double x1, double x2){
        xMin = x1;
        xMax = x2;
    }
};

class FirstOrderCDSolver{
private:
    int nNodes = 0;
    double timeStep = 0;
    double initialTime = 0;
    double finalTime = 0;
    double xStep = 0;
    pFunc mComputeCoeffA;
    pFunc mComputeCoeffC;
private:
    FirstOrderCD *mpPDE;
    BoundaryConditions_1D *mpCondition;
public:
    void SetNumberNodes(int n){
        nNodes = n;
    }
    void SetTimeStep(double t){
        timeStep = t;
    };
    void SetInitialTime(double t){
        initialTime = t;
    };
    void SetFinalTime(double t){
        finalTime = t;
    };
public:
    FirstOrderCDSolver(FirstOrderCD* pde, BoundaryConditions_1D* cond){
        mpCondition = cond;
        mpPDE = pde;
    }
public:
    bool ShouldChangeTimeStep(){
        return 0;
    }
    void CentralExplicitSolve();
    void UpwindSolve();
    void PrepareCentral();
    void PrepareUpWindSolve();
    void ComputeCentral();
    void ComputeUpWind();
    
    double UpWindFuncC(double t, double x){
        return this -> mpPDE -> mFuncC(t,x) + 1 /2.0 * mpPDE -> mFuncA(t,x) * xStep;
    }
    
};


#endif /* defined(__ConvectionDiffusion_1D__ConvectionDiffusion_1D__) */
