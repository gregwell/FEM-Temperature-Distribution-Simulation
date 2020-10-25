#pragma once
#include <cmath>

//Calculation of integral of given function f(x,y) over the interval (-1,1) using Gaussian Quadrature Method.

double f(double x, double y)
{
	return -5 * x*x*y + 2 * x*y*y + 10;
}

double inline calculate_integral_2p()
{
	double ip = 1 / sqrt(3);
	return f(-ip,-ip)+ f(ip, -ip) + f(ip, ip) + f(-ip, ip);
}

double inline calculate_integral_3p()
{
	double ip = sqrt(0.6);
	double sum = 0.0;
	double x[3];
	double w[3];

	x[0] = -ip;
	x[1] = 0.0;
	x[2] = ip;
	w[0] = 5.0/9.0;
	w[1] = 8.0/9.0;
	w[2] = w[0];
	
	for (auto i = 0; i < 3 ; i++)
	{
		for ( auto j = 0; j < 3 ; j ++)
			sum += f(x[j], x[i])*w[i] * w[j];
	}
	return sum;
}
