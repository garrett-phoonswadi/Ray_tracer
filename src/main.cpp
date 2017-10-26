//Garrett Phoonswadi
//cs410
//8/30/16
#include <model.h>
#include "imageCreator.h"
// #include <windows.h>
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;

int main(int argc, char* argv[]){
	// cout << "It's a start" << endl;
	model m;
	ifstream camera(argv[1]);
	if (camera.fail())
		return -1;
	
	m.readCameraFile(camera);
	// ifstream istr(argv[2]);
	// if (istr.fail())
	// 	return -1;
	// m.readPlyFile(istr);

	imageCreator im;
	im.createCameraVariables(&m);
	im.createCameraCorrdinateSystem(&m);
	// // im.printCameraCoord();
	im.createImagePlane();
	im.writeOut();
	// m.clearPoints();
	return 0;
}