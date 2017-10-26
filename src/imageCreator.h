#ifndef IMAGECREATOR_H_INCLUDE
#define IMAGECREATOR_H_INCLUDE

#include "model.h"
#include "matrix.h"
#include "point.h"
#include "sphere.h"
#include "light.h"
#include <cstdlib>
using std::malloc;
using std::free;
#include <algorithm>
using std::max;
#include <string>
#include <math.h>
using std::string;
#include <utility>
using std::pair;
using std::tuple;
class imageCreator{

	public: 
		imageCreator(){}

		void createCameraCorrdinateSystem(model* m);
		void createImagePlane();
		void pixelPt(double i, double j, double* p);
		void createCameraVariables(model* m);
		pair<double,int>  calculateIntersection(double* ray, double* pixpt);
		pair<double, int>  checkSphereIntersection(double* U, double* L);
		void writeOut();
		int* calculateColor(double t);
		int* ray_trace(double* ray, double* pixpt, double* accum, double* refatt, int level);
		int* calculatePolyColor(double* intersecPt, int triPos, double* L);
		double* polyHelper(int triPos);
		void printArray(double* p, int size, string s);
		void printVector(vector<double>* v);
		void printCameraCoord();
		void freeCamerCoord();
		matrix neo;
		vector< vector< pair< pair<double,int>, pair<double*,double*> > >* >* imagePlane;
		vector<vector<int*>* >* newColors;
		vector<double* > normalVectors;
		vector<point*>* points;
		vector<point*>* faces;
		vector<sphere*>* spheres;
		vector<light*>* lights;
		double width;
		double height;
		double left;
		double right;
		double top;
		double bottom;
		double near;
		double maxD;
		double minD;

		double eye [3];
		double look [3];
		double up [3];
		double ambient[3];
		double ka[3];
		double ks[3];
		double kd[3];
		double kr[3];
		double mat[3];

		double UC [3];
		double VC [3];
		double WC [3];
		string outfileName;
	private:
};

#endif //IMAGECREATOR_H_INCLUDE