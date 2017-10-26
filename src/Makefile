CXXFLAGS= -std=c++11 -O2 -Wall -g -I.

raytracer: main.o model.o point.o imageCreator.o matrix.o light.o sphere.o poly_model.o object.o
	g++ $(CXXFLAGS) main.o model.o point.o imageCreator.o matrix.o light.o sphere.o poly_model.o object.o -o raytracer
main.o: main.cpp
	g++ $(CXXFLAGS) -c main.cpp
model.o: model.cpp model.h
	g++ $(CXXFLAGS) -c model.cpp
point.o: point.cpp point.h
	g++ $(CXXFLAGS) -c point.cpp
imageCreator.o: imageCreator.cpp imageCreator.h
	g++ $(CXXFLAGS) -c imageCreator.cpp
matrix.o: matrix.cpp matrix.h
	g++ $(CXXFLAGS) -c matrix.cpp
light.o: light.cpp light.h
	g++ $(CXXFLAGS) -c light.cpp
sphere.o: sphere.cpp sphere.h
	g++ $(CXXFLAGS) -c sphere.cpp
poly_model.o: poly_model.cpp poly_model.h
	g++ $(CXXFLAGS) -c poly_model.cpp
object.o: object.cpp object.h
	g++ $(CXXFLAGS) -c object.cpp
clean:
	rm -f raytracer *.o *.ppm *.tar
	

target:
	tar cvf raytracer.tar *.cpp *.h *.txt Makefile README
