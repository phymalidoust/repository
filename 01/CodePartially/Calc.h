#pragma once
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <complex>
#include <fstream>
#include <thread>

#include "Matrix.h"
#include "DielectricResponse.h"

class Calc {
private:
	DielectricResponse DR;

	double OLow;
	double OHigh;
	double OStep;

	double wHigh;
	double wLow;
	double wStep;
	
	double thetaLow;
	double thetaHigh;
	double thetaStep;

	unsigned int OSteps;
	unsigned int wSteps;
	unsigned int thetaSteps;
	unsigned int kSteps;

	std::vector<double> xValues;
	std::vector<double> yValues;

	std::vector<std::thread> threads;
public:
	Calc(DielectricResponse DR, double OLow, double OHigh, unsigned int OSteps, unsigned int wSteps, unsigned int thetaSteps, unsigned int kSteps);

	double getwDelta() const;
	double getwLow(double O) const;
	double getwHigh(double O) const;

	double getkDelta() const;
	std::vector<double> getkInterval(double O, double w) const;

	void integrateMidpoint(double& result, double O, double i);

	double integrate(double O);
	std::vector<std::vector<double>> initializeCalc();

	//void integrateTrapez(double& result, double O, double w);
};