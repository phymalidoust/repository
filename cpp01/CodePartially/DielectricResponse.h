#pragma once
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <complex>
#include <thread>

#include "Matrix.h"

class DielectricResponse {
private:
    
//Matrix velocities
    Matrix vy;
    Matrix vx;

    double temp; //Temperature
    double kb; //Bultzman
    double u; //Potential
    double vf; //Fermi velocity
    double hBar;
    double eta; // Imaginary number

public:
    Matrix hamiltonian(double k, double theta) const;//Create the Hamiltonian
    Matrix spectral(double k, double w, double theta) const;//Create the spectral function
    double fermiFunction(double x) const;

    DielectricResponse();
    double integrand(double O, double w, double theta, double k) const;
    double getConst(char which) const;
};