//
//  ConvectionDiffusion.h
//  ConvectionDiffusion
//
//  Created by Li Xinrui on 10/25/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//


#ifndef __ConvectionDiffusion__ConvectionDiffusion__
#define __ConvectionDiffusion__ConvectionDiffusion__

typedef double (*pFunc)(double t, double x, double y);
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

#include <stdio.h>
#include "Matrix.h"

using std::cout;
using std::cin;
using std::endl;

class SecondOrderCD{
    friend class CDSolver;
    friend class CDEigenSolver;
private:
    pFunc mFuncA;
    pFunc mFuncC;
public:
    SecondOrderCD(pFunc a, pFunc c){
        mFuncA = a;
        mFuncC = c;
    };
};


class BoundaryConditions{
    friend class CDSolver;
    friend class CDEigenSolver;

private:
    bool isDirichelt = 0;
    bool isNeunmann = 0;
    pFunc mpBoundary;
    pFunc mpInitial;
    double xMin;
    double xMax;
    double yMin;
    double yMax;
public:
    void SetDirichlet(){
        isDirichelt = 1;
    }
    void SetBoundaryFunc(pFunc f){
        mpBoundary = f;
    }
    void SetInitialCondition(pFunc f){
        mpInitial = f;
    }
    void SetRegion(double x1, double x2, double y1, double y2){
        xMin = x1;
        xMax = x2;
        yMin = y1;
        yMax = y2;
    };
};

class CDSolver{
private:
    int nNodes;
    double TimeStep;
    double InitialTime;
    double FinalTime;
    double xStep;
private:
    SecondOrderCD *mpPDE;
    BoundaryConditions* mpCond;
public:
    void SetNumberNodes(int n){
        nNodes = n;
    }
    void SetTimeStep(double t){
        TimeStep = t;
    };
    void SetInitialTime(double t){
        InitialTime = t;
    };
    void SetFinalTime(double t){
        FinalTime = t;
    };
public:
    CDSolver(SecondOrderCD* pde, BoundaryConditions* cond){
        mpPDE = pde;
        mpCond = cond;
    }
public:
    bool isTimeStepLarge(){
        return 0;
    }
    void SetSpaceStep(){
        xStep = (this-> mpCond -> xMax - this -> mpCond -> xMin) / (nNodes - 1.0); // Space Step(Uniform)
    };
    
    void updateBoundary(Matrix& u, double time);
    void GenerateVector(Matrix& uPost); // Initial uPost With Values
    
    void ExplicitSolve();
    void ExplicitComputeStep(const Matrix& A, Matrix& uPre, Matrix& uPost, double AtTime);
    void GenerateExplicitMatrix(Matrix& A, double atTime); // Depending on TimeStep and nNodes;
    
    void ImplicitSolve();
    void ImplicitComputeStep(const Matrix& A, const Matrix& B, Matrix& uPre, Matrix& uPost, double atTime);
    void GenerateImplicitMatrices(Matrix& A, Matrix& B, double atTime);
    
    void SetReference(Matrix& ref, double atTime);
};

// Eigen Library, Comparison
//    void SolvePrintFile();







#endif /* defined(__ConvectionDiffusion__ConvectionDiffusion__) */
