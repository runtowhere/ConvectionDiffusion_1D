//
//  Tridiagonal.cpp
//  ConvectionDiffusion_1D
//
//  Created by Li Xinrui on 11/1/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#include <stdio.h>
#include "Matrix.h"
#include "Tridiagonal.h"

void TridiagonalSolve(mVector& a, mVector& b, mVector& c, mVector& d, mVector& sol){
    int dim = b.dim();
    int i;
    mVector cPrime(dim);
    mVector dPrime(dim);
    cPrime[0] = c[0]/b[0];
    dPrime[0] = d[0]/b[0];
    for (i = 1; i != dim - 1 ; i++) {
        assert(b[i] - a[i] * cPrime[i - 1] != 0);
        cPrime[i] = c[i] / (b[i] - a[i] * cPrime[i - 1]);
    }
    for (i = 1; i != dim; i++) {
        assert((b[i] - a[i] * cPrime[i-1]) != 0);
        dPrime[i] = (d[i] - a[i] * dPrime[i-1]) / (b[i] - a[i] * cPrime[i-1]);
    }
    sol[dim - 1] = dPrime[dim - 1];
    for(i = dim - 2; i >= 0 ; i--) {
        sol[i] = dPrime[i] - cPrime[i] * sol[i + 1];
    }
}