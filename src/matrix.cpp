#include "matrix.h"


double * matrix::subtract(double* a, double* b){
	static double s[3];
	s[0] = a[0]-b[0];
	s[1] = a[1]-b[1];
	s[2] = a[2]-b[2];

	return s;
}

double * matrix::add(double* a, double* b){
	static double s[3];
	s[0] = a[0]+b[0];
	s[1] = a[1]+b[1];
	s[2] = a[2]+b[2];

	return s;
}

double * matrix::normalize(double* a){
	static double u [3];
	double mag = sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
	u[0] = (double)a[0]/mag;
	u[1] = (double)a[1]/mag;
	u[2] = (double)a[2]/mag;

	return u;
}

double * matrix::crossProduct(double * a, double * b){
	double x = (a[1]*(double)b[2]) -(a[2]*(double)b[1]);
	double y = (a[2]*(double)b[0]) -(a[0]*(double)b[2]);
	double z = (a[0]*(double)b[1]) -(a[1]*(double)b[0]);
	static double u [3];
	u[0] = x;
	u[1] = y;
	u[2] = z;
	return u;
}

double * matrix::scale(double scaler, double* mat){
	double x = scaler*mat[0];
	double y = scaler*mat[1];
	double z = scaler*mat[2];
	static double fin [3];
	fin[0] = x;
	fin[1] = y;
	fin[2] = z;
	// cout << "x: " << x <<" y: " << y << " z:" << z << endl;
	return fin;
}

double matrix::dotProduct(double* a, double* b){
	double x = a[0]*b[0];
	double y = a[1]*b[1];
	double z = a[2]*b[2];
	return x+y+z;
}

double * matrix::pair_wise(double* a, double* b){
	static double p[3];
	p[0] = a[0]*b[0];
	p[1] = a[1]*b[1];
	p[2] = a[2]*b[2];
	return p;
}