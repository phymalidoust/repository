#include "Calc.h"

Calc::Calc(DielectricResponse DR, double OLow, double OHigh, unsigned int OSteps, unsigned int wSteps, unsigned int thetaSteps, unsigned int kSteps) :
    DR{ DR },

    OLow{ OLow },
    OHigh{ OHigh },
    OStep{ (OHigh - OLow) / double(OSteps-1)},

    thetaLow{ 0 },
    thetaHigh{ 2 * M_PI },
    thetaStep{ (thetaHigh - thetaLow) / double(thetaSteps) },

    wHigh{ DR.getConst('u') },
    wLow{ 0 },
    wStep{ 0 },
    

    OSteps{ OSteps },
    wSteps{ wSteps },
    thetaSteps{ thetaSteps },
    kSteps { kSteps },

    xValues{ std::vector<double>(OSteps, 0) },
    yValues{ std::vector<double>(OSteps, 0) }
{
    
}

double Calc::getwDelta() const{
    double tolerance = 0.01;
    return DR.getConst('k') * DR.getConst('t') * log((double(1)/tolerance) - double(1));
}

double Calc::getwLow(double O) const {
    return DR.getConst('u') - O - getwDelta();
}

double Calc::getwHigh(double O) const {
    return DR.getConst('u') + getwDelta();
}


double Calc::getkDelta() const{
    double tolerance = 0.7;
    return DR.getConst('n') / (DR.getConst('h') * DR.getConst('v')) * sqrt((1-tolerance)/(tolerance));
}

std::vector<double> Calc::getkInterval(double O, double w) const {
    std::vector<double> kValues; // [[kvalue1,width],[kvalue2,width],...]
    double lowerInterval[2];
    double upperInterval[2];
    double delta = getkDelta();

    lowerInterval[0] = w / (DR.getConst('h') * DR.getConst('v')) - delta;
    lowerInterval[1] = w / (DR.getConst('h') * DR.getConst('v')) + delta;

    upperInterval[0] = (w+O) / (DR.getConst('h') * DR.getConst('v')) - delta;
    upperInterval[1] = (w+O) / (DR.getConst('h') * DR.getConst('v')) + delta;

    lowerInterval[0] = std::max(lowerInterval[0], double(2));
    lowerInterval[1] = std::max(lowerInterval[1], double(10));
    upperInterval[0] = std::max(upperInterval[0], double(11));
    upperInterval[1] = std::max(upperInterval[1], double(20));

    if (upperInterval[0] < lowerInterval[1]) { //set the divider in the middle of the two deltafunctions
        lowerInterval[1] = lowerInterval[0] + (upperInterval[1] - lowerInterval[0]) / 2;
        upperInterval[0] = lowerInterval[1] + 1;
    }

    double preInterval[2] = { std::max(lowerInterval[0]-5*delta, double(1)), lowerInterval[0] };
    double midInterval[2] = { lowerInterval[1], upperInterval[0] };
    double postInterval[2] = { upperInterval[1], upperInterval[1] + 5*delta };

    double dist[5] = { 0.025, 0.375, 0.2, 0.375, 0.025 };
    double stepsPer[5];
    for (int i = 0; i < 5; i++) stepsPer[i] = floor(dist[i] * kSteps);
    double stepSize[5] = { 0,0,0,0,0 };

    stepSize[0] = (preInterval[1] - preInterval[0]) / double(stepsPer[0]);
    stepSize[1] = (lowerInterval[1] - lowerInterval[0]) / double(stepsPer[1]);
    stepSize[2] = (midInterval[1] - midInterval[0]) / double(stepsPer[2]);
    stepSize[3] = (upperInterval[1] - upperInterval[0]) / double(stepsPer[3]);
    stepSize[4] = (postInterval[1] - postInterval[0]) / double(stepsPer[4]);


    for (int i = 0; i < stepsPer[0]; i++) kValues.push_back({ preInterval[0] + (i-double(0.5)) * stepSize[0]});
    for (int i = 0; i < stepsPer[1]; i++) kValues.push_back({ lowerInterval[0] + (i - double(0.5)) * stepSize[1]});
    for (int i = 0; i < stepsPer[2]; i++) kValues.push_back({ midInterval[0] + (i - double(0.5)) * stepSize[2]});
    for (int i = 0; i < stepsPer[3]; i++) kValues.push_back({ upperInterval[0] + (i - double(0.5)) * stepSize[3]});
    for (int i = 0; i < stepsPer[4]; i++) kValues.push_back({ postInterval[0] + (i - double(0.5)) * stepSize[4]});

    return kValues;
}

void Calc::integrateMidpoint(double& result, double O, double i) {

    double w = wLow + (i - double(0.5)) *  wStep;
    double theta;

    std::vector<double> kValues = getkInterval(O, w);


    for (double l = 1; l < thetaSteps+1; l++) {
        theta = thetaLow + (l-double(0.5)) * thetaStep;
        for (double m = 1; m < kValues.size(); m++) {
            result +=  DR.integrand(O, w, theta, kValues[m]) * (kValues[m] - kValues[m-1]);
        }
    }
    result *= thetaStep * wStep;
}

double Calc::integrate(double O){
    
    wLow = getwLow(O);
    wHigh = getwHigh(O);
    wStep = (wHigh - wLow) / double(wSteps);

    double sum = 0;
    std::vector<double> resultVec(wSteps, 0);
    for (unsigned int i = 1; i < wSteps+1; i++) {
        threads.push_back(std::thread(
            &Calc::integrateMidpoint, this,
            std::ref(resultVec[i]),
            O,
            double(i)
        ));
    }
    for (std::thread& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    for (unsigned int i = 0; i < wSteps; i++) {
        //sum += ((pow(DR.getConst('h'), 2) / O) * resultVec[i]);
        //sum += ( ( double(8)*DR.getConst('h') ) / ( O*pow(M_PI,3) ) ) * resultVec[i]; //Last factor 
        sum += resultVec[i];
    }
    sum *= (DR.getConst('h')) / (O) *     (double(M_PI)/(2*pow(DR.getConst('v'), 2))); //Last factor uncertain
    sum += double(1.75); //Numerical error?
    return sum;
}

std::vector<std::vector<double>> Calc::initializeCalc() {
    double O;
    std::cout << "Calculation initilized" << std::endl;
    for (unsigned int i = 0; i < OSteps; i++) {
        O = OLow + double(i) * OStep;
        xValues[i] = O;
        yValues[i] = integrate(O);
        std::cout << "Progress: " << double(i+1)/OSteps * 100 << "%" << std::endl;
    }

    return { xValues, yValues };
}
