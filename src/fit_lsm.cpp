/*
 * fit_lsm.cpp
 *
 *  Created on: 22 июн. 2017 г.
 *      Author: user
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <cmath>
#include "fit_lsm.h"

void linearFitLeastSquare(const std::vector<double>& y, const std::function<double(int)> x, double& a, double& b){

	int n = y.size();
	double xsum=0, x2sum=0, ysum=0, xysum=0;
    for (int i = 0; i <  n; i++)
    {
        xsum = xsum+x(i);
        ysum = ysum+y[i];
        x2sum = x2sum+pow(x(i),2);
        xysum = xysum+x(i)*y[i];
    }
    a = (n * xysum - xsum * ysum)/(n * x2sum - xsum * xsum);            //calculate slope
    b = (x2sum * ysum - xsum * xysum)/(x2sum * n - xsum * xsum);            //calculate intercept
}

void polinomialFitLeastSquare(const std::vector<double>& y, const std::function<double(int)> x, const int degree, std::vector<double>& poli){

	int N = y.size();
	int n = degree;
	    std::vector<double> X;                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	    for (int i = 0; i < 2 * n + 1; i++)
	    {
	    	X.push_back(0);
	        for (int j = 0; j < N; j++){
	            X.back() += pow(x(j), i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	        }
	    }


	    std::vector<std::vector<double> > B;
	    for (int i = 0; i <= n; i++){
	    	B.push_back(std::vector<double>());
	    	for (int j = 0; j <= n; j++){
	    		B.back().push_back(X[i+j]);
	    	}
	    }

	    std::vector<double> Y;                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	    for (int i = 0; i < n + 1; i++)
	    {
	        Y.push_back(0);
	        for (int j = 0; j < N; j++){
	        	Y.back()+= pow(x(j), i) * y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	        }
	    }

	    for (int i = 0; i <= n; i++){
	        B[i][n + 1] = Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
	    }

	    n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

//	    std::cout<<"\nThe Normal(Augmented Matrix) is as follows:\n";
//	    for (int i=0;i<n;i++)            //print the Normal-augmented matrix
//	    {
//	        for (int j=0;j<=n;j++)
//	        	std::cout<<B[i][j]<< " ";
//	        std::cout<<"\n";
//	    }

	    for (int i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
	        for (int k=i+1;k<n;k++)
	            if (B[i][i]<B[k][i])
	                for (int j=0;j<=n;j++)
	                {
	                    double temp=B[i][j];
	                    B[i][j]=B[k][j];
	                    B[k][j]=temp;
	                }

	    for (int i=0;i<n-1;i++)            //loop to perform the gauss elimination
	        for (int k=i+1;k<n;k++)
	            {
	                double t=B[k][i]/B[i][i];
	                for (int j=0;j<=n;j++)
	                    B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
	            }


	    std::vector<double> a(n, 0);
	    for (int i=n-1;i>=0;i--)                //back-substitution
	    {                        //x is an array whose values correspond to the values of x,y,z..
	        a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
	        for (int j=0;j<n;j++)
	            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
	                a[i]=a[i]-B[i][j]*a[j];
	        a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	    }
	    poli = std::move(a);
}



