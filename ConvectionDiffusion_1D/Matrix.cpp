//
//  Matrix.cpp
//  MatrixClass
//
//  Created by Li Xinrui on 10/1/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "Matrix.h"
#include <cmath>
using std::cout;
using std::endl;
using std::cin;

// Algorithms
int LUDecomp(const Matrix& m, Matrix& l, Matrix& u){ // l, u blank Matrices
    int n = m.row();
    if(m.row() != m.col()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    int i, j, k; // Copy
    for(int i = 0; i != n; ++i)
        for(int j = 0; j != n; ++j){
            u[i][j] = m[i][j];
        }
    
    for(i = 0; i != n; ++i) l[i][i] = 1;
    
    for(k = 0; k != n - 1; ++k){ // L
        for(i = k + 1; i != n; ++i){
            l[i][k] = u[i][k] / u[k][k];
        }
        for(i = k + 1; i != n; ++i){ // U
            u[i][k] = 0;
            for(j = k + 1; j != n; ++j){
                u[i][j] = u[i][j] - l[i][k] * u[k][j];
            }
        }
    }
    return 0;
}

int LUPivotDecomp(const Matrix& m, Matrix& p, Matrix& l, Matrix& u){ // Mathematical Subscripts

    int n = m.row();
    int i, j, k;
    int *per = new int[n - 1]; // at k, exchange row k with row per[k - 1]
    TYPE max = 0;
    TYPE temp = 0;
    //Test
    if(m.row() != m.col()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    //Copy M into U
    for(i = 0; i != n; ++i)
        for(j = 0; j != n; ++j){
            u[i][j] = m[i][j];
        }
    //Initiate L
    for(i = 0; i != n; ++i) l[i][i] = 1;

    //Computation
    for(k = 1; k <= n; k++){
        // Find Permutations
        max = std::abs(u[k-1][k-1]);
        per[k-1] = k;
        for(i = k; i <= n; i++){
            if(std::abs(u[i-1][k-1]) > max){
                max = std::abs(u[i-1][k-1]);
                per[k-1] = i;
            }
        }
        u.rowSwap(k,per[k - 1]);
        // Compute L
        for(i = k + 1; i <= n; i++){
            l[i-1][k-1] = u[i-1][k-1] / u[k-1][k-1];
        }
        // Update U
        for(i = k + 1; i <= n; i++){
            u[i - 1][k - 1] = 0;
            for(j = k + 1; j <= n; ++j){
                u[i - 1][j - 1] = u[i - 1][j - 1] - l[i - 1][k - 1] * u[k - 1][j - 1];
            }
        }
    }
    // Compute P
    for(i = 1; i <= n; i++) p[i-1][i-1] = 1;
    for(k = 1; k <= n-1; k++) p.rowSwap(k, per[k-1]);
    
    // Compute L
    for(k = 1; k != n - 1; k++){
        if(per[k] != k + 1 ){
            for(j = 1; j != k + 1; j++){  // Exchange Row k+1 and row per[k] in the Lower Left Matrix
                temp = l[k][j-1];
                l[k][j-1] = l[ (per[k] - 1) ][j-1];
                l[ (per[k] - 1) ][j-1] = temp;
            }
        }
    }

    return 0;
}

int CholeskyDecomp(const Matrix& m, Matrix& l){
    int i, p, k;
    int n = m.row();
    TYPE temp;
    
    l[0][0] = sqrt(m[0][0]);
    for(i = 1; i <= n; i++){
        l[i - 1][0] = m[i - 1][0] / l[0][0];
    }
    for(k = 2; k <= n; k++){
        temp = 0;
        for(p = 1; p <= k - 1; p++){
            temp += l[k - 1][p - 1] * l[k - 1][p - 1];
        }
        if(m[k - 1][k - 1] - temp < 0){ // Negative Diagonal 
            cout << "Cholesky: Diagonal cannot Sqrt at " << k << endl;
            return 1;
        }
        else{
        l[k - 1][k - 1] = sqrt(m[k - 1][k - 1] - temp);
        }
        for(i = k + 1; i <= n; i++){
            temp = 0;
            for(p = 1; p <= k - 1; p++){
                temp += l[i - 1][p - 1] * l[k - 1][p - 1];
            }
            l[i - 1][k - 1] = (m[i - 1][k - 1] - temp) / l[k-1][k-1];
        }
    }
    
    return 0;
}

int LDLDecomp(const Matrix& m, Matrix& l, Matrix& d){
    int i, j, k;
    int n = m.row();
    TYPE temp;
    for(i = 1; i <= n; i++){
        l[i - 1][i - 1] = 1;
    }
    for(j = 1; j <= n; j++){
        TYPE* v = new TYPE[j - 1];
        temp = 0;
        for(k = 1; k <= j - 1; k++){
            v[k - 1] = d[k - 1][k - 1] * l[j - 1][k - 1];
            temp += l[j - 1][k - 1] * v[k - 1];
        }
        d[j - 1][j - 1] = m[j - 1][j - 1] - temp;
        for(i = j + 1; i <= n; i++){
            temp = 0;
            for(k = 1; k <= j - 1; k++){
                temp += l[i -1][k - 1] * v[k - 1];

            }
            l[i - 1][j - 1] = (m[i-1][j-1] - temp) / d[j-1][j-1];
        }
        delete v;
    }
    return 0;
}


int LUSolve(const Matrix& m, const Matrix& rhs, Matrix& res){ // NEED RECYCLE
    if(rhs.col() == 1){
        int n = m.row();
        Matrix l(n,n), u(n,n);
        LUDecomp(m, l, u);
        Matrix temp(n,1);
        lowerSolve(l, rhs, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.row();
        Matrix l(n,n), u(n,n);
        LUDecomp(m, l, u);
        for(int i = 1; i <= rhs.col(); i++){
            Matrix temp1(n,1);
            Matrix temp2(n,1);
            Matrix right(n,1);
            for(int j = 1; j <= rhs.row(); j++){
                right[j - 1][0] = rhs[j - 1][i - 1];
            }
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.row(); j++){
                res[j - 1][i - 1] = temp2[j - 1][0];
            }
        }
    }
    return 0;
}

int LUPivotSolve(const Matrix& m, const Matrix& rhs, Matrix& res){ // NEED
    if(rhs.col() == 1){
        int n = m.row();
        Matrix p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        Matrix newb(n,1);
        newb = p * rhs;
        Matrix temp(n,1);
        lowerSolve(l, newb, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.row();
        Matrix p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        for(int i = 1; i <= rhs.col(); i++){
            Matrix right(n,1);
            for(int j = 1; j <= rhs.row(); j++){
                right[j-1][0] = rhs[j - 1][i - 1];
            }
            right = p * right;
            Matrix temp1(n,1);
            Matrix temp2(n,1);
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.row(); j++){
                res[j - 1][i - 1] = temp2[j - 1][0];
            }
        }
    }
    return 0;
}

int CholeskySolve(const Matrix&m, const Matrix& rhs, Matrix& res){
    int n = m.row();
    Matrix y(n,1);
    Matrix l(n,n);
    Matrix trans(n,n);
    if(CholeskyDecomp(m, l)){
        cout << "Cholesky Solve: Cannot Decomp" << endl;
        return 1;
    }
    lowerSolve(l, rhs, y);
    l.transpose(trans);
    upperSolve(trans, y, res);
    return 0;
}
int LDLSolve(const Matrix& m, const Matrix& rhs, Matrix& res){
    int n = m.row();
    Matrix y(n,1);
    Matrix l(n,n), d(n,n);
    Matrix trans(n,n);
    LDLDecomp(m, l, d);
    l.transpose(trans);
    lowerSolve(l, rhs, y);
    upperSolve(d * trans, y, res);
    return 0;
}

int lowerSolve(const Matrix& l,const Matrix& rhs, Matrix& res) // Array or Vector? I love vector.. Without Check on Dim
{
    int n = l.row();
    res[0][0] = rhs[0][0] / l[0][0];
    for(int i = 1; i != n; ++i){
        TYPE temp = 0;
        for(int j = 0; j <= i - 1; ++j){
            temp += l[i][j] * res[j][0];
        }
        res[i][0] = (rhs[i][0] - temp) / l[i][i];
    }
    return 0;
}

int upperSolve(const Matrix& u, const Matrix& rhs, Matrix& res){
    int n = u.row();
    res[n - 1][0] = rhs[n - 1][0] / u[n - 1][n - 1];
    for(int i = n - 2; i >= 0; --i){
        TYPE temp = 0;
        for(int j = i + 1; j != n; ++j){
            temp += u[i][j] * res[j][0];
        }
        res[i][0] = (rhs[i][0] - temp) / u[i][i];
    }
    return 0;
}

// Destructors
Matrix::~Matrix(){
    int i;
    for (i = 0; i != nRow; i++) {
        delete [] mat[i];
    }
    delete [] mat;
}

// Constructors
Matrix::Matrix():
nRow(0), nCol(0), mat(nullptr) {
}

Matrix::Matrix(int row, int col):
nRow(row), nCol(col), mat(new TYPE*[row]) {
    for(int i = 0; i != row; ++i){
        mat[i] = new TYPE[col];
    }
    for(int i = 0; i != row; ++i)
        for(int j = 0; j!= col; ++j){
            mat[i][j] = 0;
        }
}

Matrix::Matrix(int dim):
nRow(dim), nCol(dim), mat(new TYPE*[dim]) {
    for(int i = 0; i != dim; ++i){
        mat[i] = new TYPE[dim];
    }
    for(int i = 0; i != dim; ++i)
        for(int j = 0; j!= dim; ++j){
            if(i == j) mat[i][j] = 1;
            else mat[i][j] = 0;
        }
}

Matrix::Matrix(const Matrix& m):
nRow(m.row()),nCol(m.col()),mat(new TYPE*[m.row()]){
    for(int i = 0; i != m.row(); ++i){
        mat[i] = new TYPE[m.col()];
    }
    for(int i = 0; i != m.row(); ++i)
        for(int j = 0; j!= m.col(); ++j){
            mat[i][j] = m[i][j];
        }
}

// Destructor

// Operators

TYPE* Matrix::operator[](int index){
    return mat[index];
}

const TYPE* Matrix::operator[](int index) const{
    return mat[index];
}

TYPE Matrix::operator()(int rowIndex, int colIndex){
    return mat[rowIndex][colIndex];
}



Matrix& Matrix::operator+=(const Matrix& rhs){
    if(this->nCol != rhs.col() || this->nRow != rhs.row()){
        std::cout << "Unequal Dimensions" << std::endl;
    }
    for(size_t i = 0; i != rhs.row(); ++i){
        for(size_t j = 0; j != rhs.col(); ++j){
            mat[i][j] += rhs.mat[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix& rhs){
    return *this += rhs * -1;
}

Matrix& Matrix::operator*=(const TYPE& rhs){
    for(int i = 0; i != this->nRow; i++)
        for(int j = 0; j != nCol; j++) {
            mat[i][j] *= rhs;
        }
    return *this;
}

Matrix& Matrix::operator/=(const TYPE& rhs){
    if(rhs == 0) std::cout << "Divide by Zero" << std::endl;
    return *this *= 1/rhs ;
}

Matrix& Matrix::operator=(const TYPE* rhs){
    int k = 0;
    for(int i = 0; i != this->nRow; i++)
        for(int j = 0; j != this->nCol; j++) {
            mat[i][j] = rhs[k];
            k++;
        }
    return *this;
}

Matrix& Matrix::operator=(const vector<TYPE>& rhs){
    int k = 0;
    for(int i = 0; i != this->nRow; i++)
        for(int j = 0; j != this->nCol; j++) {
            mat[i][j] = rhs[k];
            k++;
        }
    return *this;
}
Matrix& Matrix::operator=(const Matrix& rhs){
    int i,j;
    if(mat == nullptr){
        this->nRow = rhs.row();
        this->nCol = rhs.col();
        mat = new TYPE*[rhs.row()];
        for(int i = 0; i != rhs.row(); ++i){
            mat[i] = new TYPE[rhs.col()];
        }
    }
    if(this->nCol != rhs.col() || this->nRow != rhs.row()){
        cout << "unequal dimensions cannot assign matrix" << endl;
    }
    else{
        for(i = 0; i != rhs.row(); i++)
            for(j = 0; j != rhs.col(); j++){
                mat[i][j] = rhs[i][j];
            }
    }
    return *this;
}


Matrix& Matrix::rowSwap(int i, int j){
    if(i != j){
    TYPE *p = mat[i-1];
    mat[i-1] = mat[j-1];
    mat[j-1] = p;
    }
    return *this;
}

Matrix& Matrix::colSwap(int i, int j) { // Mathematical Subscripts
    TYPE temp;
    for(int k = 0; k != this->nRow; k++){
        temp = mat[k][i - 1];
        mat[k][i - 1] = mat[k][j - 1];
        mat[k][j - 1] = temp;
    }
    return *this;
}
bool Matrix::isSquare(){
    if(this->nRow == this->nCol){
        return 1;
    }
    else return 0;
}

double Matrix::Norm2Vec(){
    int i,j;
    double sum = 0;
    for (i = 0; i != this->nRow; i++) {
        for (j = 0; j != this->nCol; j++) {
            sum += mat[i][j] * mat[i][j];
        }
    }
    return sqrt(sum);
}

double Matrix::NormInfVec(){
    double max = 0;
    for (int i = 0; i != this->nRow; i++) {
        if (std::abs((double) this->mat[i][0]) > max ) {
            max = std::abs((double) this->mat[i][0]);
        }
    }
    return max;
};

Matrix& Matrix::transpose(){
    int i, j;
    TYPE temp;
    // This changes the Matrix itself
    if(this->nRow != this->nCol){
        Matrix copy(*this);
        delete mat;
        mat = new TYPE*[this->nCol];
        for(i = 0; i != this->nCol; i++){
            mat[i] = new TYPE[this->nRow];
        }
        for(i = 1; i <= nRow; i++)
            for(j = 1; j <= nCol; j++){
                mat[j - 1][i - 1] = copy[i - 1][j - 1];
            }
        this->nCol = copy.row();
        this->nRow = copy.col();
        return *this;
    }
    for(i = 1; i <= this->nRow; i++)
        for(j = 1; j < i; j++){
            if(mat[i - 1][j - 1] != mat[j - 1][i - 1]){
                temp = mat[i - 1][j - 1];
                mat[i - 1][j - 1] = mat[j - 1][i - 1];
                mat[j - 1][i - 1] = temp;
            }
        }
    return *this;
}

int Matrix::transpose(Matrix& res){
    int i, j;
    for(i = 1; i <= this->nRow; i++)
        for(j = 1; j <= this->nCol; j++){
            res[j - 1][i - 1] = mat[i - 1][j - 1];
        }
    return 0;
}


// Public

Matrix operator+(const Matrix& m1, const Matrix& m2){
    Matrix sum(m1.row(), m1.col());
    sum = m1;
    sum += m2;
    return sum;
}

Matrix operator-(const Matrix& m1, const Matrix& m2){
    Matrix sum(m1.row(), m1.col());
    sum = m1;
    sum -= m2;
    return sum;
}
Matrix operator*(const Matrix& m1, const Matrix& m2){
    if(m1.col() != m2.row()) std::cout << "Illegal Multiplication" << std::endl;
    Matrix res(m1.row(), m2.col());
    for(int i = 0; i != m1.row(); ++i)
        for(int j = 0; j != m2.col(); ++j)
            for(int k = 0; k != m1.col(); ++k){
                res[i][j] += (m1[i][k] * m2[k][j]);
            }
    return res;
}
Matrix operator*(const Matrix& m1, const TYPE& num){
    Matrix res(m1.row(),m1.col());
    res = m1;
    res *= num;
    return res;
}
Matrix operator/(const Matrix& m1, const TYPE& num){
    Matrix res(m1.row(),m1.col());
    res = m1;
    res /= num;
    return res;
}


std::ostream& operator<<(std::ostream& out, const Matrix& m){
    int i, j;
    if(m.mat == nullptr) cout << "Empty Matrix";
    for (i = 0; i != m.row(); i++){
        for(j = 0; j!= m.col(); j++){
            out << m[i][j] << " ";
        }
        out << std::endl;
    }
    return out;
}

std::ostream& latex(std::ostream& out, const Matrix& m, int line){
    int i, j;
    int t = 0;
    if(m.mat == nullptr) cout << "Empty Matrix";
    for (i = 0; i != m.row(); i++){
        for(j = 0; j!= m.col(); j++){
            if(i == m.row() - 1 && j == m.col() - 1){
                out << m[i][j] << " " << "\\\\";
            }
            else{
                t++;
                if(t == line){
                    out << m[i][j] << " ";
                    out << "\\\\";
                    out << endl;
                    out << "\\hline";
                    out << endl;
                    t = 0;
                }
                else{
                    out << m[i][j] << " " << "&" << " ";
                }
            }
        }
    }
    return out;
}



