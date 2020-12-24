#pragma once
#include <iostream>

using namespace std;

class Gauss {
private:

	int n; //the number of equations
	double **a;
	double *x;

public:
	void gaussElimination(int n_n);

	~Gauss() {
		//delete[] x, y;
	}
};