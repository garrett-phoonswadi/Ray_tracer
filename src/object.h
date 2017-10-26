#ifndef OBJECT_H_INCLUDE
#define OBJECT_H_INCLUDE


class object{
	double x,y,z,w;
	public: object(double x, double y, double z, double w){

	}

	inline void setX(double a){x=a;}
	inline double getX(){return x;}

	inline void setY(double b){y =b;}
	inline double getY(){return y;}

	inline void setZ(double c){z=c;}
	inline double getZ(){return z;}

	inline void setW(double d){w =d;}
	inline double getW(){return w;}

};
#endif
