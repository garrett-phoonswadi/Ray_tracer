//Garrett Phoonswad
// create a model from a specific directions
#include <model.h>

void model::readPlyFile(istream& istr){
	string reader;
	string temp;
	bool found_endheader = false;
	bool found_points = false;
	int vertex = 0;
	// int face =0;
	int vertices =0;
	string rounded = "_rounded";
	string centered = "_centered";
	while(getline(istr, reader)){
		// cout << "input string " <<reader << endl;
		stringstream input_string(reader);
		//Create header vector
		while(input_string >> temp && !found_endheader){
			// cout << "Start of string " << temp << endl;
			if (temp.compare("element") == 0){
				input_string >> temp;
				// cout << "middle of string " << temp << endl;
				if (temp.compare("vertex") == 0){
					input_string >>temp;
					// cout << "vetex found " << temp <<endl;
					vertex = atoi(temp.c_str());
				}
				if (temp.compare("face") ==0){
					input_string >>temp;
					// cout <<"face found " <<temp <<endl;
					face = atoi(temp.c_str());
				}
			}
			// cout << "input string " << reader << endl;

		}
		if (reader.compare("end_header") ==0){
			found_endheader = true;
			header.push_back(reader);
			continue;
			}
		if (!found_endheader){
			header.push_back(reader);

		}
		if (found_endheader)
			vertices++;
		if (found_endheader && vertices <= vertex){
			// new code
			int start = 0;
			int next = 0;
			string tempx;
			string tempy;
			string tempz;
			while (reader.at(0) == ' ')
				reader = reader.substr(1, reader.size());
			next = reader.find(" ") +1;
			tempx = reader.substr(start, next);
			reader = reader.substr(next, reader.size());
			// cout << "x " << tempx << endl;
			// cout << "readerY:" << reader << endl;
			while (reader.at(0) == ' ')
				reader = reader.substr(1, reader.size());
			// cout << "readerY:" << reader << endl;
			next = reader.find(" ");
			tempy = reader.substr(start, next);
			// cout << "y " << tempy << endl;
			reader = reader.substr(next, reader.size());
			// cout << "readerZ:" << reader << endl;
			while (reader.at(0) == ' ')
				reader = reader.substr(1, reader.size());
			// reader = reader.substr(next, reader.size());
			next = reader.size();
			tempz = reader.substr(start, next);
			// cout << "z " << tempz << endl;
			// cout << "x " << tempx << " y " << tempy << " z " << tempz<< endl;
			// create point object from strings
			double x = 0.0;
			double y = 0.0;
			double z = 0.0;
			x = atof(tempx.c_str());
			y = atof(tempy.c_str());
			z = atof(tempz.c_str());
			point* p = new point(x,y,z);
			points.push_back(p);
		}
		if (vertices == vertex){
			found_points = true;
			continue;
		}
		if (found_points && vertices > vertex){
			int start = 0;
			int next = 0;
			string numSide;
			string tempx;
			string tempy;
			string tempz;
			next = reader.find(" ") +1;
			numSide = reader.substr(start, next);
			reader = reader.substr(next, reader.size());

			next = reader.find(" ") +1;
			tempx = reader.substr(start, next);
			reader = reader.substr(next, reader.size());
			// cout << "x " << tempx << endl;
			next = reader.find(" ") +1;
			tempy = reader.substr(start, next);
			// cout << "y " << tempy << endl;
			reader = reader.substr(next, reader.size());
			next = reader.size();
			tempz = reader.substr(start, next);
			// cout << "z " << tempz << endl;
			// cout << "x " << tempx << " y " << tempy << " z " << tempz<< endl;
			// create point object from strings
			double x = 0.0;
			double y = 0.0;
			double z = 0.0;
			x = atof(tempx.c_str());
			y = atof(tempy.c_str());
			z = atof(tempz.c_str());
			// cout << "a: " << x << " b: " << y << " c: " << z << endl;
			point* f = new point(x,y,z);
			faces.push_back(f);

		}

	}
	// cout << "points---------------" << endl;
	// printPoints(points);
	// cout << "faces----------------" << endl;
	// printPoints(faces);

}

void model::printPoints(vector<point*>& v){
	// cout << "printing points" << endl;

	for (unsigned int i =0; i < v.size(); i++){
		double x = v[i]->getX();
		double y = v[i]->getY();
		double z = v[i]->getZ();
		cout <<  x  << " " << y << " " << z << endl;
	}
}

void model::readCameraFile(istream& istr){

	string reader;
	string temp;
	while(getline(istr, reader)){
		stringstream input_string(reader);
		while(input_string >> temp){
			if (temp == "eye"){
				input_string >> temp;
				double x = atof(temp.c_str());
				eye[0] = x;
				input_string >> temp;
				double y = atof(temp.c_str());
				eye[1] = y;
				input_string >> temp;
				double z = atof(temp.c_str());
				eye[2] = z;
			}
			else if (temp == "look"){
				input_string >> temp;
				double x = atof(temp.c_str());
				look[0] = x;
				input_string >> temp;
				double y = atof(temp.c_str());
				look[1] = y;
				input_string >> temp;
				double z = atof(temp.c_str());
				look[2] = z;
			}
			else if (temp == "up"){
				input_string >> temp;
				double x = atof(temp.c_str());
				up[0] = x;
				input_string >> temp;
				double y = atof(temp.c_str());
				up[1] = y;
				input_string >> temp;
				double z = atof(temp.c_str());
				up[2] = z;

			}
			else if (temp == "d"){
				input_string >> temp;
				distance = atof(temp.c_str());
			}
			else if (temp == "bounds"){
				input_string >> temp;
				double x = atof(temp.c_str());
				bounds[0] = x;
				input_string >> temp;
				double y = atof(temp.c_str());
				bounds[1] = y;
				input_string >> temp;
				double z = atof(temp.c_str());
				bounds[2] = z;
				input_string >> temp;
				double zi = atof(temp.c_str());
				bounds[3] = zi;
			}
			else if (temp == "res"){
				input_string >> temp;
				double x = atof(temp.c_str());
				res[0] = x;
				input_string >> temp;
				double y = atof(temp.c_str());
				res[1] = y;

			}
			else if (temp == "ambient"){
				input_string >> temp;
				ambient[0] = atof(temp.c_str());
				input_string >> temp;
				ambient[1] = atof(temp.c_str());
				input_string >> temp;
				ambient[2] = atof(temp.c_str());
			}
			else if (temp == "light"){
				// cout << "reading in lgihts" << endl;
				input_string >> temp;
				double x = atof(temp.c_str());
				input_string >> temp;
				double y = atof(temp.c_str());
				input_string >> temp;
				double z = atof(temp.c_str());
				input_string >> temp;
				double w = atof(temp.c_str());
				input_string >> temp;
				double r = atof(temp.c_str());
				input_string >> temp;
				double g = atof(temp.c_str());
				input_string >> temp;
				double b = atof(temp.c_str());
				// cout << "x: " << x << " y: " << y << " z: " << z << endl;
				light* l = new light(x,y,z,w,r,g,b);
				lights.push_back(l);

			}
			else if (temp == "sphere"){
				input_string >> temp;
				double x = atof(temp.c_str());
				input_string >> temp;
				double y = atof(temp.c_str());
				input_string >> temp;
				double z = atof(temp.c_str());
				input_string >> temp;
				double radius = atof(temp.c_str());
				input_string >> temp;
				double r = atof(temp.c_str());
				input_string >> temp;
				double g = atof(temp.c_str());
				input_string >> temp;
				double b = atof(temp.c_str());
				sphere* s = new sphere(x,y,z,radius,r,g,b);
				spheres.push_back(s);

			}
			else if (temp == "model"){
				input_string >> temp;
				double x = atof(temp.c_str());
				input_string >> temp;
				double y = atof(temp.c_str());
				input_string >> temp;
				double z = atof(temp.c_str());
				input_string >> temp;
				// axis angle 
				double ax = atof(temp.c_str()); 
				input_string >> temp;
				double ay = atof(temp.c_str());
				input_string >> temp;
				double az = atof(temp.c_str());
				input_string >> temp;
				double at = atof(temp.c_str());
				input_string >> temp;
				poly_model* p = new poly_model(x,y,z,ax,ay,az,at,temp);
				poly.push_back(p);
			}
			else if (temp == "ks"){
					input_string >> temp;
					ks[0] = atof(temp.c_str());
					input_string >> temp;
					ks[1] = atof(temp.c_str());
					input_string >> temp;
					ks[2] = atof(temp.c_str());
				}
			else if (temp == "ka"){
				input_string >> temp;
				ka[0] = atof(temp.c_str());
				input_string >> temp;
				ka[1] = atof(temp.c_str());
				input_string >> temp;
				ka[2] = atof(temp.c_str());
			}
			else if (temp == "kd"){
				input_string >> temp;
				kd[0] = atof(temp.c_str());
				input_string >> temp;
				kd[1] = atof(temp.c_str());
				input_string >> temp;
				kd[2] = atof(temp.c_str());
			}
			else if (temp == "kr"){
				input_string >> temp;
				kr[0] = atof(temp.c_str());
				input_string >> temp;
				kr[1] = atof(temp.c_str());
				input_string >> temp;
				kr[2] = atof(temp.c_str());
			}
			else if (temp == "outfile"){
				input_string >> temp;
				outfile = temp;
			}
	}

	}
	// cout << "camera file" << endl;
	// double* p = eye;
	// printArray(p, 3);
	// p = look;
	// printArray(p, 3);
	// p = up;
	// printArray(p, 3);
	// cout << distance << endl;
	// p = bounds;
	// printArray(p, 4);
	// p = res;
	// printArray(p, 2);
	// printLights();
	// printPoly();
	// cout << "end camera file" << endl;
	readObjFiles();
}

void model::readObjFiles(){
	string temp;
	string line;
	for(unsigned int i=0; i< poly.size(); i++){
		string filename = poly[i]->getFileName();
		cout << "fileName: " << filename << endl;
		ifstream istr(filename);
		while (getline(istr,temp)){
			stringstream input_string(temp);
			while(input_string >> line){
				if (line == "v"){ // points
					// cout << "input: " << temp << endl;
					input_string >> line;
					double x = atof(line.c_str());
					input_string >> line;
					double y = atof(line.c_str());
					input_string >> line;
					double z = atof(line.c_str());
					cout << "x: " << x << " y: " << y << " z: " << z << endl;
					point* p = new point(x,y,z);
					points.push_back(p);
				}
				else if(line == "f"){// faces
					input_string >> line;
					double x = atof(line.c_str());
					input_string >> line;
					double y = atof(line.c_str());
					input_string >> line;
					double z = atof(line.c_str());
					point* f = new point(x,y,z);
					faces.push_back(f);
				}
				else if (line == "ks"){
					input_string >> line;
					ks[0] = atof(line.c_str());
					input_string >> line;
					ks[1] = atof(line.c_str());
					input_string >> line;
					ks[2] = atof(line.c_str());
				}
				else if (line == "ka"){
					input_string >> line;
					ka[0] = atof(line.c_str());
					input_string >> line;
					ka[1] = atof(line.c_str());
					input_string >> line;
					ka[2] = atof(line.c_str());
				}
				else if (line == "kd"){
					input_string >> line;
					kd[0] = atof(line.c_str());
					input_string >> line;
					kd[1] = atof(line.c_str());
					input_string >> line;
					kd[2] = atof(line.c_str());
				}
				else if (line == "mat"){
					input_string >> line;
					mat[0] = atof(line.c_str());
					input_string >> line;
					mat[1] = atof(line.c_str());
					input_string >> line;
					mat[2] = atof(line.c_str());
				}
				else{
					break;
				}
			}
		}
	}
	// printObj();
}

void model::printArray(double* p, int size){
	for (int i = 0; i < size; i++){
		cout << p[i] << " ";
	}
	cout << endl;
}

void model::clearPoints(){
	for (unsigned int i =0; i < points.size(); i++){
		delete points.at(i);
	}

	for (unsigned int i =0; i < faces.size(); i++){
		delete faces.at(i);
	}
}

void model::printLights(){
	cout << "printing lights" << endl;
	for (unsigned int i =0; i<lights.size(); i++){
		light* l = lights[i];
		cout << l->getX() << " " << l->getY() << " " << l->getZ() << endl;
		cout << l->getW() << endl;
		cout << l->getR() << " " << l->getG() << " " << l->getB() << endl;
		cout << "---------" << endl;
	}
}

void model::printPoly(){
	cout << "print poly" << endl;
	for (unsigned int i =0; i< poly.size(); i++){
		poly_model* p = poly[i];
		cout << "wx: " << p->getWorldX() << " wy: " << p->getWorldY() << " wz: " << p->getWorldZ() << endl;
		cout << "ax: " << p->getAxisX() << " ay: " << p->getAxisY() << " az: " << p->getAxisZ() << " at: " << p->getAxisT() << endl;
		cout << "filename: " << p->getFileName() << endl;
		cout << "------------------------------" << endl;
	}
}

void model::printObj(){
	for (unsigned int i=0; i< objFiles.size(); i++){
		vector<object*> temp = objFiles[i];
		for(unsigned int j=0; j<temp.size(); i++){
			cout << "x: " << temp[j]->getX() << " y: " << temp[j]->getY() << " z: " << temp[j]->getZ() << endl;
		}
	}
}