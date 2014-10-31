//
//  ConvectionDiffusion.cpp
//  ConvectionDiffusion
//
//  Created by Li Xinrui on 10/25/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "ConvectionDiffusion.h"
#include "Matrix.h"
#include <vector>
using std::cout;
using std::cin;
using std::endl;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;


void CDSolver::updateBoundary(Matrix& u, double time){
    int i,j;
    int node = this -> nNodes;
    for (i = 0; i != node; i++) {
        u[i * node + 0][0] = mpCond -> mpBoundary(time, 0 * xStep, i * xStep);
        u[i * node + (node-1)][0] = mpCond -> mpBoundary(time, (node-1) * xStep, i * xStep);
    }
    for (j = 0; j != node; j++) {
        u[0 * node + j][0] = mpCond -> mpBoundary(time, j * xStep, 0 * xStep);
        u[(node-1) * node + j][0] = mpCond -> mpBoundary(time, j * xStep, (node-1) * xStep);
    }
}

void CDSolver::ExplicitComputeStep(const Matrix& A, Matrix& uPre, Matrix& uPost, double AtTime){
    uPre = uPost;
    updateBoundary(uPre, AtTime);
    uPost = A * uPre;
}

void CDSolver::ExplicitSolve(){
    SetSpaceStep();
    Matrix A(nNodes * nNodes, nNodes * nNodes);
    Matrix uPre(nNodes * nNodes,1);
    Matrix uPost(nNodes * nNodes,1);
    GenerateVector(uPost);
    int i;
    double tNow;
    int shouldBreak = 0;
    for(i = 1; (tNow = (i * TimeStep + InitialTime)) <= FinalTime; i++){
        if(FinalTime - tNow < TimeStep && (tNow - FinalTime) > 0.001 * TimeStep){
            TimeStep = FinalTime - tNow;
            shouldBreak = 1;
        }
        GenerateExplicitMatrix(A,tNow);
        ExplicitComputeStep(A, uPre, uPost, tNow);
        if (isTimeStepLarge()) {
//            SetTimeStep(compute from time and space);
        }
        if (shouldBreak) {
            break;
        }
    }
    // Reference
    Matrix ref(this -> nNodes * nNodes,1);
    SetReference(ref, FinalTime);
    cout << (ref - uPost).NormInfVec() << ",";
//    cout << (ref - uPost).NormInfVec();
}

void CDSolver::GenerateExplicitMatrix(Matrix& A, double atTime){
    int i, j;
    int node = this -> nNodes;
    for(i = 0; i != node; i++){
        for(j = 0; j != node; j++){
            double u = this->mpPDE->mFuncC(atTime, j * xStep, i * xStep) * TimeStep / (xStep * xStep);
            double v = this->mpPDE->mFuncA(atTime, j * xStep, i * xStep) * TimeStep / (xStep * 2);
            if(i == 0 || i == node - 1 || j == 0 || j == node - 1){ // On Boundary
                A[i * node + j][i * node + j] = 1; //  updating Boundary!!!!
            }
            else{
                A[i * node + j][(i+1) * node + j] = u - v;
                A[i * node + j][i * node + (j+1)] = u - v;
                A[i * node + j][i * node + j] = 1 - 4 * u;
                A[i * node + j][i * node + (j-1)] = u + v;
                A[i * node + j][(i-1) * node + j] = u + v;
            }
        }
    }
}

void CDSolver::GenerateVector(Matrix& uPost){
    int i,j;
    double node = this -> nNodes;
    for(i = 0; i != node; i++){
        for(j = 0; j != node; j++){
            uPost[i * node + j][0] = this->mpCond->mpInitial(InitialTime, j * xStep, i * xStep);
        }
    }
}

// ImplicitSolver

void CDSolver::GenerateImplicitMatrices(Matrix& A, Matrix& B, double atTime){
    int i, j;
    int node = this->nNodes;
    for(i = 0; i != node; i++){
        for(j = 0; j != node; j++){
            double u = this->mpPDE->mFuncC(atTime, j * xStep, i * xStep) * TimeStep / (2 * xStep * xStep);
            double v = this->mpPDE->mFuncA(atTime, j * xStep, i * xStep) * TimeStep / (4 * xStep);
            if(i == 0 || i == node - 1 || j == 0 || j == node - 1){ // On Boundary
                A[i * node + j][i * node + j] = 1; //  updating Boundary!!!!
                B[i * node + j][i * node + j] = 1;
            }
            else{
                // AAAAA (v-u, v-u, 1 + 4u, -(u+v),-(v+u))
                A[i * node + j][(i+1) * node + j] = v - u;
                A[i * node + j][i * node + (j+1)] = v - u;
                A[i * node + j][i * node + j] = 1 + 4 * u;
                A[i * node + j][i * node + (j-1)] = -(u + v);
                A[i * node + j][(i-1) * node + j] = -(u + v);
                // BBBBB (u-v,u-v, 1-4u, u+v, u+v)
                B[i * node + j][(i+1) * node + j] = u - v;
                B[i * node + j][i * node + (j+1)] = u - v;
                B[i * node + j][i * node + j] = 1 - 4 * u;
                B[i * node + j][i * node + (j-1)] = (u + v);
                B[i * node + j][(i-1) * node + j] = (u + v);
            }
        }
    }
//    cout << A;
}

void CDSolver::ImplicitComputeStep(const Matrix& A, const Matrix& B, Matrix& uPre, Matrix& uPost, double atTime){
    uPre = uPost;
    updateBoundary(uPre, atTime);
    LUPivotSolve(A, B * uPre, uPost);
}


void CDSolver::ImplicitSolve(){
    SetSpaceStep();
    Matrix A(nNodes * nNodes, nNodes * nNodes);
    Matrix B(nNodes * nNodes, nNodes * nNodes);
    Matrix uPre(nNodes * nNodes,1);
    Matrix uPost(nNodes * nNodes,1);
    GenerateVector(uPost);
    int i;
    int shouldBreak = 0;
    double tNow = 0;
    GenerateImplicitMatrices(A, B, tNow);
    for (i = 1; (tNow = i * TimeStep + InitialTime) <= FinalTime; i++) {
        if(FinalTime - tNow < TimeStep && (tNow - FinalTime) > 0.001 * TimeStep){
            TimeStep = FinalTime - tNow;
            shouldBreak = 1;
        }
        GenerateImplicitMatrices(A, B, tNow);
        ImplicitComputeStep(A, B, uPre, uPost, tNow);
        if (isTimeStepLarge()) {
//            SetTimeStep(newTimeStep);
        }
        if (shouldBreak) {
            break;
        }
    }
    Matrix ref(this->nNodes * this->nNodes,1);
    SetReference(ref, FinalTime);
    cout << (ref - uPost).NormInfVec() << ",";
//    cout << (ref - uPost).NormInfVec();
}

void CDSolver::SetReference(Matrix& ref, double atTime){
    int i, j;
    for(i = 0; i != this -> nNodes; i++){
        for (j = 0; j != this -> nNodes; j++) {
            ref[i * nNodes + j][0] = this -> mpCond -> mpBoundary(atTime, j * xStep, i * xStep);
        }
    }
}
// EIGENSOLVER USING SPARSE ARRAYS

void CDEigenSolver::GenerateVector(Eigen::VectorXd& u){
    int i,j;
    double node = this -> nNodes;
    for(i = 0; i != node; i++){
        for(j = 0; j != node; j++){
            u(i * node + j) = this->mpCond->mpInitial(InitialTime,mpCond->xMin + j * xStep, mpCond->yMin + i * xStep);
        }
    }
}

void CDEigenSolver::updateBoundary(Eigen::VectorXd& u, double atTime){
    int i,j;
    int node = this -> nNodes;
    for (i = 0; i != node; i++) {
        u(i * node + 0) = mpCond -> mpBoundary(atTime, 0 * xStep, i * xStep);
        u(i * node + (node-1)) = mpCond -> mpBoundary(atTime, (node-1) * xStep, i * xStep);
    }
    for (j = 0; j != node; j++) {
        u(0 * node + j) = mpCond -> mpBoundary(atTime, j * xStep, 0 * xStep);
        u((node-1) * node + j) = mpCond -> mpBoundary(atTime, j * xStep, (node-1) * xStep);
    }
}
void CDEigenSolver::EigenComputeStep(const SpMat& A, const SpMat& B, Eigen::VectorXd& uPre, Eigen::VectorXd& uPost, double atTime,
                                     Eigen::SparseLU<Eigen::SparseMatrix<double> >& solver){
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        cout << "Decomp Fail" << endl;
        return;
    }
    uPre = uPost;
    updateBoundary(uPre, atTime);
    uPost = solver.solve(B * uPre);
    if (solver.info() != Eigen::Success) {
        cout << "Solve Fail" << endl;
    }
}

void CDEigenSolver::EigenSolve(){
    SetSpaceStep();
    Eigen::VectorXd uPre(nNodes * nNodes);
    Eigen::VectorXd uPost(nNodes * nNodes);
    GenerateVector(uPost);
    SpMat A(nNodes * nNodes, nNodes * nNodes);
    SpMat B(nNodes * nNodes, nNodes * nNodes);
    int i;
    double tNow;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    for (i = 1; (tNow = i * TimeStep + InitialTime) <= FinalTime; i++) {
        GenerateMatrix(A, B, tNow);
        EigenComputeStep(A, B, uPre, uPost, tNow, solver);
    }
    Eigen::VectorXd ref(nNodes * nNodes);
    SetReference(ref, FinalTime);
    ref -= uPost;
    cout << EigenVecNormInf(ref) << endl;
}

void CDEigenSolver::GenerateMatrix(SpMat &A, SpMat &B, double atTime){
    std::vector<T> aTripleList;
    std::vector<T> bTripleList;
    aTripleList.reserve(nNodes * nNodes);
    bTripleList.reserve(nNodes * nNodes);
    int i, j;
    for (i = 0; i != nNodes; i++) {
        for (j = 0; j != nNodes; j++) {
            double u = this->mpPDE->mFuncC(atTime, j * xStep, i * xStep) * TimeStep / (2 * xStep * xStep);
            double v = this->mpPDE->mFuncA(atTime, j * xStep, i * xStep) * TimeStep / (4 * xStep);
            if (i == 0 || i == nNodes - 1 || j == 0 || j == nNodes - 1) { // Boundary
                aTripleList.push_back(T(i * nNodes + j, i * nNodes + j, 1.0));
                bTripleList.push_back(T(i * nNodes + j, i * nNodes + j, 1.0));
            }
            else{
                // AAAAA
                aTripleList.push_back(T(i * nNodes + j, (i+1) * nNodes + j, v - u));
                aTripleList.push_back(T(i * nNodes + j, i * nNodes + (j+1), v - u));
                aTripleList.push_back(T(i * nNodes + j, i * nNodes + j, 1 + 4 * u));
                aTripleList.push_back(T(i * nNodes + j, (i) * nNodes + (j-1), -(u+v)));
                aTripleList.push_back(T(i * nNodes + j, (i-1) * nNodes + (j), -(u+v)));
                // BBBBB
                bTripleList.push_back(T(i * nNodes + j, (i+1) * nNodes + j, u - v));
                bTripleList.push_back(T(i * nNodes + j, i * nNodes + (j+1), u - v));
                bTripleList.push_back(T(i * nNodes + j, i * nNodes + j, 1 - 4 * u));
                bTripleList.push_back(T(i * nNodes + j, (i) * nNodes + (j-1), (u+v)));
                bTripleList.push_back(T(i * nNodes + j, (i-1) * nNodes + (j), (u+v)));
            }
        }
    }
    A.setFromTriplets(aTripleList.begin(), aTripleList.end());
    B.setFromTriplets(bTripleList.begin(), bTripleList.end());
}

void CDEigenSolver::SetReference(Eigen::VectorXd &ref, double atTime){
    int i, j;
    for(i = 0; i != this -> nNodes; i++){
        for (j = 0; j != this -> nNodes; j++) {
            ref(i * nNodes + j) = this -> mpCond -> mpBoundary(atTime, j * xStep, i * xStep);
        }
    }
}
