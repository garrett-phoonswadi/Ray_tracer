#ifndef MODEL_H_INCLUDE
#define MODEL_H_INCLUDE

// includes
#include "point.h"

#include <vector>
using std::vector;
#include <string>
using std::string;
// using std::stoi;
// using std::substr;
#include <fstream>
using std::istream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
using std::stringstream;
#include <sstream>
using std::istringstream;
#include <iostream>
using std::cout;
using std::endl;
using std::ifstream;
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <cstring>
using std::strcpy;
#include <array>
#include "light.h"
#include "poly_model.h"
#include "sphere.h"
#include "object.h"

class model{
	public: model(){}

		void readPlyFile(istream& istr);
		void printPoints(vector<point*>& v);
		void printArray(double* p, int size);
		void calcualteStandardDivation(vector<point*>& p, double* mv);
		void calcualteMean(vector<point*>& p, double* sd);
		void getBoundingBox(vector<point*>& p);
		void centerPoints(vector<point*>& p, double* mv);
		void whitenPoints(vector<point*>& p, double* sd);
		void printLights();
		void printPoly();
		void printObj();
		void readCameraFile(istream& istr);
		void readObjFiles();
		void clearPoints();
		inline vector<point*>* getPoints() { return &points;}
		inline vector<point*>* getFaces() { return &faces;}
		inline vector<sphere*>* getSpheres() {return &spheres;}
		inline vector<light*>* getLights(){return &lights;}

		inline double* getEye() { return eye;}
		inline double* getLook() { return look;}
		inline double* getUp() { return up;}
		inline double* getBounds() { return bounds;}
		inline double* getRes() { return res;}
		inline double getDistance() { return distance;}
		inline double* getAmbient(){return ambient;}
		inline double* getKa(){return ka;}
		inline double* getKd(){return kd;}
		inline double* getKs(){return ks;}
		inline double* getKr(){return kr;}
		inline double* getMat(){return mat;}
		inline string getOutfileName(){return outfile;}
		//  vectors for ply file
		vector<string> header;
		vector<point*> points;
		vector<point*> faces;
		vector<vector<object*> > objFiles;
		vector<light*> lights;
		vector<poly_model*> poly;
		vector<sphere*> spheres;
		int face;
		// vectors for camera
		double eye[3];
		double look[3];
		double up[3];
		double distance;
		double bounds [4];
		double res[2];
		double ambient[3];
		double ks[3];
		double ka[3];
		double kd[3];
		double kr[3];
		double mat[3];
		string outfile;
		

	private:
		

};

#endif // MODEL_H_INCLUDE