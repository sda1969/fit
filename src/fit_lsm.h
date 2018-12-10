/*
 * fit_lsm.h
 *
 *  Created on: 22 июн. 2017 г.
 *      Author: user
 */

#ifndef FIT_LSM_H_
#define FIT_LSM_H_

#include <vector>
#include <functional>

/* Линейная и Полиномиальная аппроксимация методом наименьших квадратов:
 * задается вектор значений y и функция вычисляющая значение х по порядковому номеру пары (x,y)
 * например лямбда [](int i){return i;}
 * linearFitLeastSquare вычисляет а - коэффициент наклона b - смещение
 * polinomialFitLeastSquare - необходимо задать дополнительно степень полинома degree
 * вычисляется вектор poli - коэффициенты полинома начиная с младшей степени
 *
 * контрольный пример:
 *
  y {1 1.8 1.3 2.5 6.3}
  для polinomialFitLeastSquare должно получиться:
   The values of the coefficients are as follows:
     x^0=1.42
     x^1=-1.07
     x^2=0.55
  для linear FitLeastSquare должно получиться:
   The values of the coefficients are as follows:
	x^0=0.32
	x^1=1.13

 * http://www.bragitoff.com/2015/09/c-program-for-polynomial-fit-least-squares/
 *
 *	polinomialFitLeastSquare можно использовать и для линейной аппроксимации, указав degree=1
 *	но будет работать медленнее чем linearFitLeastSquare
 *	7701 ns time interval for linear
 *	8757 ns time interval for polinomial
 * */

void linearFitLeastSquare(		const std::vector<double>& y,
								const std::function<double(int)> x,
								double& a, double& b);

void polinomialFitLeastSquare(	const std::vector<double>& y,
								const std::function<double(int)> x,
								const int degree,
								std::vector<double>& poli);




#endif /* FIT_LSM_H_ */
