//============================================================================
// Name        : front.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <cmath>
#include <functional>
#include <chrono>

#include "fit_lsm.h"
#include <unistd.h>

using namespace std::chrono;

int main(int argc, char * argv[]) {
	std::vector<double> points;

	//typedef duration < long long unsigned int, std::ratio<1,1000000000> > llnanoseconds;

	printf("number of points = %d\n", argc - 1);
	for(int i = 1; i < argc; ++i){
		points.push_back(atof(argv[i]));
		if (points.back() == 0){
			printf("argument[%d] is invalid!\n", i);
			return -1;
		}
		printf("point[%d] = %f\n", i - 1, points.back());
	}

	double a, b;
	//steady_clock::time_point t1 = steady_clock::now();
	auto t1 = steady_clock::now();
    linearFitLeastSquare(points, [](int i){return  i;}, a, b);
    auto t2 = steady_clock::now();
    nanoseconds  dur = duration_cast < nanoseconds > (t2 - t1);
    std::cout << dur.count() << " ns time interval for linear \n";
    printf("y = %5.3fx + %5.3f\n",  a, b);

    std::vector<double> poli;
    auto t3 = steady_clock::now();
    polinomialFitLeastSquare(points, [](int i){return i;}, 2, poli);
    auto t4 = steady_clock::now();
    nanoseconds  dur2 = duration_cast < nanoseconds >(t4 - t3);
    std::cout << dur2.count() << " ns time interval for polinomial\n";

    std::cout<<"\nThe values of the coefficients are as follows:\n";
    for (size_t i = 0; i < poli.size(); i++){
    	std::cout<<"x^"<<i<<"="<<poli[i]<<std::endl;            // Print the values of x^0,x^1,x^2,x^3,....
    }
    std::cout<<"\n";



	return 0;
}
