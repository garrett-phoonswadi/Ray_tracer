#ifndef POINT_H_INCLUDE
#define POINT_H_INCLUDE

class point{
	double x, y, z;
	public: 

	point(double a, double b, double c){
		x=a;
		y=b;
		z=c;
	}

	inline double getX() const {
		return x;
	}

	inline double getY() const{
		return y;
	}

	inline double getZ() const{
		return z;
	}
	inline void setX(double a) {
		x=a;
	}

	inline void setY(double a){
		y=a;
	}

	inline void setZ(double a){
		z=a;
	}


	private:
};

#endif //POINT_H_INCLUDE