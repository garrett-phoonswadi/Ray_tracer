#ifndef LIGHT_H_INCLUDE
#define LIGHT_H_INCLUDE

#include <iostream>
using std::cout;
using std::endl;


class light{
	double x,y,z,w,r,g,b;
	public: light(double a, double m, double c, double d, double e, double f, double l){
		// cout << "x: " << a << " y: " << m << " z: " << c << endl;
		// cout << "w: " << d << endl;
		// cout << "r: " << e << " g: " << f << " b: " << l << endl;
		setX(a);
		setY(m);
		setZ(c);
		setW(d);
		setR(e);
		setG(f);
		setB(l);
		// cout << "sx: " << getX() << " sy: " << getY() << " sz: " << getZ() << endl;
	}

		inline double getX() const { return x;}
		inline double getY() const { return y;}
		inline double getZ() const { return z;}
		inline double getW() const { return w;}
		inline double getR() const { return r;}
		inline double getG() const { return g;}
		inline double getB() const { return b;}

		inline void setX(double a) {x = a;}
		inline void setY(double m) {y = m;}
		inline void setZ(double c) {z = c;}
		inline void setW(double d) {w = d;}
		inline void setR(double e) {r = e;}
		inline void setG(double f) {g = f;}
		inline void setB(double l) {b = l;}

		inline double* getPos(){
			static double pos [3];
			pos[0] = getX();
			pos[1] = getY();
			pos[2] = getZ();
			return pos;
		}

		inline double* getColor(){
			static double color [3];
			color[0] = getR();
			color[1] = getG();
			color[2] = getB();
			return color;
		}

	private:


};
#endif // MODEL_H_INCLUDE