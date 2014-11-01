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
#include "Tridiagonal.h"

using std::cout;
using std::cin;
using std::endl;

typedef double (*pFunc)(double t, double x, double u);
typedef double (*pFunc_2)(double t, double x);
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
    pFunc_2 mpBoundaryFunc;
    pFunc_2 mpInitialFunc;
    double xMin;
    double xMax;
public:
    BoundaryConditions_1D(pFunc_2 boundary, pFunc_2 initial){
        mpBoundaryFunc = boundary;
        mpInitialFunc = initial;
    }
public:
    void SetDirichlet(){
        isDirichelt = 1;
    }
    void SetBoundaryFunc(pFunc_2 f){
        mpBoundaryFunc = f;
    }
    void SetInitialFunc(pFunc_2 f){
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
    void SetXStep(){
        xStep = (mpCondition -> xMax - mpCondition -> xMin) / (nNodes - 1);
    }
    void UpdateBoundary(mVector& u, double atTime);
    void SetInitialValue(mVector& u);
    
    void CentralExplicitSolve();
    void ComputeCentral(mVector& uPre, mVector& uPost, double atTime);
    
    void UpwindSolve();
    void ComputeUpWind(mVector& uPre, mVector& uPost, double atTime);
   
// Error when trying to define local function and pFunc to it , How to Fix?
//    double UpWindFuncC(double t, double x){
//        return this -> mpPDE -> mFuncC(t,x) + 1 /2.0 * mpPDE -> mFuncA(t,x) * xStep;
//    }
    
};





#endif /* defined(__ConvectionDiffusion_1D__ConvectionDiffusion_1D__) */
