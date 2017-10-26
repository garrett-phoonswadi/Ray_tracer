#ifndef MATRIX_H_INCLUDE
#define MATRIX_H_INCLUDE
#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

class matrix{

	public: 
		matrix(){}

		double * subtract(double* a, double* b);
		double * add(double* a, double* b);
		double * normalize(double* a);
		double * crossProduct(double * a, double * b);
		double * scale(double scaler, double* mat);
		double   dotProduct(double* a, double* b);
		double * pair_wise(double* a, double* b);
	private:
};

#endif //MATRIX_H_INCLUDE