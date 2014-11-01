//
//  main.cpp
//  ConvectionDiffusion_1D
//
//  Created by ikon on 10/31/14.
//  Copyright (c) 2014 ___Imaginaire___. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
#include "ConvectionDiffusion_1D.h"
using std::cout;
using std::endl;
using std::cin;

double u(double t, double x){
    return 5 * t + x;
}
pFunc p1 = u;

double f(double t, double x){
    return p1(t,x) + 1;
}


int main(int argc, const char * argv[])
{
    pFunc p1;
    p1 = u;
    pFunc p2;
    cout << f(3,4);
    
    cout << "Je suis ton pere" << endl;
    return 0;
}

