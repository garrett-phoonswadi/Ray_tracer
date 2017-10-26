#ifndef POLY_MODEL_H_INCLUDE
#define POLY_MODLE_H_INCLUDE

#include <string>
using std::string;

class poly_model{
	double wx, wy, wz, ax, ay, az, at;
	string filename;
	public: poly_model(double wx, double wy, double wz, 
		double ax, double ay, double az, double at, string file){
		setWorldX(wx);
		setWorldY(wy);
		setWorldZ(wz);
		setAxisX(ax);
		setAxisY(ay);
		setAxisZ(az);
		setAxisT(at);
		setFileName(file);

	}

	inline void setWorldX(double x){wx = x;}
	inline void setWorldY(double y){wy = y;}
	inline void setWorldZ(double z){wz = z;}
	inline void setAxisX(double a){ ax = a;}
	inline void setAxisY(double b){ ay = b;}
	inline void setAxisZ(double c){ az = c;}
	inline void setAxisT(double t){ at = t;}
	inline void setFileName(string s){ filename = s;}

	inline double getWorldX(){return wx;}
	inline double getWorldY(){return wy;}
	inline double getWorldZ(){return wz;}
	inline double getAxisX(){ return ax;}
	inline double getAxisY(){ return ay;}
	inline double getAxisZ(){ return az;}
	inline double getAxisT(){ return at;}
	inline string getFileName(){ return filename;}

};
#endif 