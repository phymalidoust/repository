#define _USE_MATH_DEFINES
//#include <catch2/catch.hpp>
#include <stdio.h> // include the input/output related functions
#include <iostream> //  responsible for input/output stream
#include <fstream> // Stream class to write on files
#include <complex> // real and imaginary numbers

#include "Matrix.h" // Containing the matrix construction operations
#include "DielectricResponse.h" // Containing the dielectric definitions
#include "Calc.h" // Performing the main calculations

int main() {
    double OLow = 0.01; // lower omega frequency
    double OHigh = 1.5; // Captital omega frequency
    unsigned int OSteps = 100;

    unsigned int wSteps = 300;
    unsigned int thetaSteps = 10;
    unsigned int kSteps = 500;

    std::cout << "Using Parameters:" << std::endl;
    std::cout << "Steps for O: " << OSteps << ", running from " << OLow << " to " << OHigh << std::endl;
    std::cout << "Steps for w:" << wSteps << std::endl;
    std::cout << "Steps for theta: " << thetaSteps << std::endl;
    std::cout << "Steps for k: " << kSteps << std::endl;

    DielectricResponse DR;
    Calc calculation(DR, OLow, OHigh, OSteps, wSteps, thetaSteps, kSteps);

    std::vector<std::vector<double>> result;
    result = calculation.initializeCalc();

    std::ofstream outputFile;
    outputFile.open("results.txt");
    for (unsigned int i = 0; i < OSteps; i++) {
        std::cout << result[0][i] << "\t" << result[1][i] << std::endl;
        outputFile << result[0][i] << "\t" << result[1][i] << std::endl;
    }
    outputFile.close();

    std::cout << "Calculation finished. Press a key to exit." << std::endl;
    std::cin.get();


    return 0;
}