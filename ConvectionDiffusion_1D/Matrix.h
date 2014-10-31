//
//  Matrix.h
//  MatrixClass
//
//  Created by Li Xinrui on 10/1/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __MatrixClass__Matrix__
#define __MatrixClass__Matrix__
#define TYPE double
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>

using std::vector;

class Matrix{
    //  Members
private:
    int nRow;
    int nCol;
    TYPE **mat;
public:
    int row(){
        return nRow;
    }
    int col(){
        return nCol;
    }
    int col() const{
        return nCol;
    }
    int row() const{
        return nRow;
    }
public:
    
    // Operators
    TYPE* operator[](int index);
    const TYPE* operator[](int index) const;
    TYPE operator()(int rowIndex, int colIndex);

    // Operations
    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    Matrix& operator*=(const TYPE& rhs);
    Matrix& operator/=(const TYPE& rhs);
    Matrix& operator=(const TYPE* rhs);
    Matrix& operator=(const vector<TYPE>& rhs);
    Matrix& operator=(const Matrix& rhs);


public:
    
    Matrix& rowSwap(int i, int j);
    Matrix& colSwap(int i, int j);
    Matrix& transpose();
    double Norm2Vec();
    double NormInfVec();
    int transpose(Matrix& res);
    bool isSquare();

    // Constructors
    Matrix();
    Matrix(int row, int col);
    Matrix(int dim);
    Matrix(const Matrix& m);
    // Deconstructor
    ~Matrix();
    
public:
    friend std::ostream& operator<<(std::ostream& out, const Matrix& m);
    friend std::ostream& latex(std::ostream& out, const Matrix& m, int line);

};

// Operations
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const TYPE& num);
Matrix operator/(const Matrix& m1, const TYPE& num);



// Algorithms

int LUDecomp(const Matrix& m, Matrix& l, Matrix& u);
int LUPivotDecomp(const Matrix& m, Matrix& p, Matrix& l, Matrix& u);
int CholeskyDecomp(const Matrix& m, Matrix& l);
int LDLDecomp(const Matrix& m, Matrix &l, Matrix& d);

int LUSolve(const Matrix& m,  const Matrix& rhs, Matrix& res);
int LUPivotSolve(const Matrix& m, const Matrix& rhs, Matrix& res);
int CholeskySolve(const Matrix&m, const Matrix& rhs, Matrix& res);
int LDLSolve(const Matrix& m, const Matrix& rhs, Matrix& res);

int lowerSolve(const Matrix& l,const Matrix& rhs, Matrix& res);
int upperSolve(const Matrix& u,const Matrix& rhs, Matrix& res);


#endif /* defined(__MatrixClass__Matrix__) */
