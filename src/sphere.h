#ifndef SPHERE_H_INCLUDE
#define SPHERE_H_INCLUDE


class sphere{
	double x,y,z,radius,r,g,b;
	public: sphere(double x, double y, double z, double rad, double red, double green, double blue){
		setX(x);
		setY(y);
		setZ(z);
		setRadius(rad);
		setR(red);
		setG(green);
		setB(blue);

	}
	inline double* getCenter(){
		static double center [3];
		center[0] = getX();
		center[1] = getY();
		center[2] = getZ();
		return center;
	}
	inline double* getColor(){
		static double color [3];
		color[0] = getR();
		color[1] = getG();
		color[2] = getB();
		return color;
	}
	inline void setX(double a){ x=a;}
	inline void setY(double b){ y=b;}
	inline void setZ(double c){ z=c;}
	inline void setRadius(double rad){ radius=rad;}
	inline void setR(double red){ r=red;}
	inline void setG(double green){ g=green;}
	inline void setB(double blue){ b=blue;}

	inline double getX(){ return x;}
	inline double getY(){ return y;}
	inline double getZ(){ return z;}
	inline double getRadius(){ return radius;}
	inline double getR(){ return r;}
	inline double getG(){ return g;}
	inline double getB(){ return b;}
};
#endif