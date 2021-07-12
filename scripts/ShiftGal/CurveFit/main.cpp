#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include <string>
#include <cmath>
#include "string.h"

using namespace alglib;

int NUM_PARAMS = 5;

void gaussian(const real_1d_array& c, const real_1d_array& x, double& func, void *ptr) 
{ 
	/*
	each value in x represents a point [x, y]
	c is a list of parameters
    c[0] = Amplitude
    c[1], c[2] = x0, y0
    c[3] = gamma
    c[4] = alpha
    */
    func = c[0] * pow(1.0 + ((pow(x[0] - c[1], 2) + pow(x[1] - c[2], 2)) / pow(c[3], 2)), -1.0 * c[4]);
    //func = c[4] * exp(-1 * pow(((x[0] - c[0]) / c[2]), 2) - pow(((x[1] - c[1]) / c[3]), 2));
}


int main(int argc, char* argv[])
{
    
	if (argc != 6) return -1;
    //if (argc != 4) return -1;

	try
    {
        real_2d_array x = argv[1];
	    real_1d_array y = argv[2];
	    real_1d_array c = argv[3];

        // bounds for each parameter
        real_1d_array boundLower = argv[4];
        real_1d_array boundUpper= argv[5];
        
        //real_1d_array boundUpper = "[10.0, 10.0, 100.0, 100.0, 10000.0]";
        //real_1d_array boundLower = "[0.0, 0.0, 0.0, 0.0, 0.0]";
        
        double epsx = 0.000001;
        ae_int_t maxits = 0;
        ae_int_t info;
        lsfitstate state;
        double diffstep = 0.0001;
        lsfitreport rep;
    
        lsfitcreatef(x, y, c, diffstep, state);
        lsfitsetbc(state, boundLower, boundUpper);
        lsfitsetcond(state, epsx, maxits);
        alglib::lsfitfit(state, gaussian);
        lsfitresults(state, info, c, rep);
   
        // output is read by python script
        std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " " << c[4]  << std::endl;
        //std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << std::endl;
    
    }
    catch (...)
    {
        return -1;
    }
	return 0;
}
