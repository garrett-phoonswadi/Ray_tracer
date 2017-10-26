#include "imageCreator.h"


void imageCreator::createCameraVariables(model* m){
	double* eyeP = m->getEye();
	eye[0] = eyeP[0];
	eye[1] = eyeP[1];
	eye[2] = eyeP[2];
	double* lookP = m->getLook();
	look[0] = lookP[0];
	look[1] = lookP[1];
	look[2] = lookP[2];
	double* upP = m->getUp();
	up[0] = upP[0];
	up[1] = upP[1];
	up[2] = upP[2];
	width = m->getRes()[0];
	height = m->getRes()[1];
	left = m->getBounds()[0];
	right = m->getBounds()[3];
	bottom = m->getBounds()[0];
	top = m->getBounds()[2];
	near = m->getDistance();
	points = m->getPoints();
	faces = m->getFaces();
	maxD = 0;
	minD = 0;
	spheres = m->getSpheres();
	lights = m->getLights();
	double* am = m->getAmbient();
	ambient[0] = am[0];
	ambient[1] = am[1];
	ambient[2] = am[2];
	double* kss = m->getKs();
	ks[0] = kss[0];
	ks[1] = kss[1];
	ks[2] = kss[2];
	double* kaa = m->getKa();
	ka[0] = kaa[0];
	ka[1] = kaa[1];
	ka[2] = kaa[2];
	double* kdd = m->getKd();
	kd[0] = kdd[0];
	kd[1] = kdd[1];
	kd[2] = kdd[2];
	double* krr = m->getKr();
	kr[0] = krr[0];
	kr[1] = krr[1];
	kr[2] = krr[2];
	double* matt = m->getMat();
	mat[0] = matt[0];
	mat[1] = matt[1];
	mat[2] = matt[2];
	newColors = new vector<vector<int*>*>();
	outfileName = m->getOutfileName();
}

void imageCreator::createCameraCorrdinateSystem(model* m){
	// cout << "in createCameraCorrdinateSystem ----------------" << endl;
	double* const eyeP = eye;
	double* const lookP = look;
	double* cameraCorrd = neo.subtract(eyeP,lookP);
	// cout << cameraCorrd[0] << " ," << cameraCorrd[1] << " ," << cameraCorrd[2] << endl;
	double* W = neo.normalize(cameraCorrd);
	WC[0] = W[0];
	WC[1] = W[1];
	WC[2] = W[2];
	// cout << W[0] << " ," << W[1] << " ," << W[2] << endl;
	double* const wp = WC;
	double* U = neo.crossProduct(up, wp);
	// cout << U[0] << " ," << U[1] << " ," << U[2] << endl;
	U = neo.normalize(U);
	UC[0] = U[0];
	UC[1] = U[1];
	UC[2] = U[2];
	// cout << U[0] << " ," << U[1] << " ," << U[2] << endl;
	double* const upt = UC;
	double* V = neo.crossProduct(wp,upt);
	VC[0] = V[0];
	VC[1] = V[1];
	VC[2] = V[2];
	// cout << VC[0] << " ," << VC[1] << " ," << VC[2] << endl;
}

void imageCreator::createImagePlane(){
	// cout << "creating image plane ---------------" << endl;
	// cout << "width: " << width << " height: " << height <<endl;
	double* const ept = eye;
	for (int j=0; j<height;j++){
		vector<int*>* temp = new vector<int*>();
		temp->reserve(width);
		for(int i=0; i<width;i++){
			double* pixpt = new double[3];
			pixelPt(i,j,pixpt);
			// cout<< "------------------pixpt: " << pixpt[0] << " " << pixpt[1] << " " << pixpt[2] << endl;
			// calulate ray
			double* const pixp = pixpt;
			double* ray = neo.subtract(pixp,ept);
			ray = neo.normalize(ray); // direction
			double* color = new double[3];
			color[0] = 0.0;
			color[1] = 0.0;
			color[2] = 0.0;
			double* refatt = new double[3];
			refatt[0] = 1.0;
			refatt[1] = 1.0;
			refatt[2] = 1.0;
			double* direction = new double[3];
			direction[0] = ray[0];
			direction[1] = ray[1];
			direction[2] = ray[2];
			double* dpt = direction;
			// cout << "x: " << i << " y: " << j << endl;
			int* coco = ray_trace(dpt, pixp, color, refatt, 3);
			int* tColor = new int[3];
			tColor[0] = coco[0];
			tColor[1] = coco[1];
			tColor[2] = coco[2];
			temp->push_back(tColor);
			delete color;
			delete refatt;
			delete direction;
		}
		newColors->push_back(temp);
	}
}

void imageCreator::pixelPt(double i, double j, double* pixpt){
	
	// cout << "i: " << i << " j: " << j << endl;
	// cout << "right: " << right << " left: " << left <<endl;
	double px = i/(width-1)*(right-left)+left;
	double py = j/(height-1)*(top-bottom)+bottom;
	// cout << "px: " << px << " py: " << py << endl;
	double* nWV;
	double nn [3];
	double* pxUV;
	double pxn[3];
	double* pvVV;
	double pyn[3];
	double* const wpt = WC;
	double* const upt = UC;
	double* const vpt = VC;
	// near*W
	nWV = neo.scale(-near, wpt);
	nn[0] = nWV[0];
	nn[1] = nWV[1];
	nn[2] = nWV[2];
	double* const bs = nn;
	// px*U
	pxUV = neo.scale(px, upt);
	pxn[0] = pxUV[0];
	pxn[1] = pxUV[1];
	pxn[2] = pxUV[2];
	double* const pxp=pxn;
	// py*V
	pvVV = neo.scale(py, vpt);
	pyn[0] = pvVV[0];
	pyn[1] = pvVV[1];
	pyn[2] = pvVV[2];
	double* const pyp=pyn;
	// printArray(bs,3);
	// printArray(pxp,3);
	// printArray(pyp,3);
	// U+V = UV
	// cout << "Adding -------------------" << endl;
	pxUV = neo.add(pxp, pyp);
	pxn[0] = pxUV[0];
	pxn[1] = pxUV[1];
	pxn[2] = pxUV[2];
	// printArray(pxp,3);
	//E+W = EW
	nWV = neo.add(bs,eye);
	nn[0] = nWV[0];
	nn[1] = nWV[1];
	nn[2] = nWV[2];
	// printArray(bs,3);
	// EW+UV = pixpt
	double* pix = neo.add(bs,pxp);
	// printArray(pix,3);
	pixpt[0] = pix[0];
	pixpt[1] = pix[1];
	pixpt[2] = pix[2];
	// delete nWV;
	// delete pxUV;
	// delete pvVV;

}

pair<double,int> imageCreator::calculateIntersection(double* d, double* pixpt){
	// cout << "calculate intersection ------------------------------------" << endl;
	double* const ept = pixpt;
	double* const dp = d;
	vector<double> tv;
	tv.reserve(faces->size());
	double min = 1000;
	int minFace =-1;
	// cout << tv.size()<<endl;
	for (unsigned int run=0; run< faces->size(); run++){
		point* const triangle = faces->at(run);
		double pointA = triangle->getX();
		double pointB = triangle->getY();
		double pointC = triangle->getZ();

		point* const pA = points->at(pointA);
		double A [3];
		A[0] = pA->getX();
		A[1] = pA->getY();
		A[2] = pA->getZ();
		point* const pB = points->at(pointB);
		double B [3];
		B[0] = pB->getX();
		B[1] = pB->getY();
		B[2] = pB->getZ();
		point* const pC = points->at(pointC);
		double C[3];
		C[0] = pC->getX();
		C[1] = pC->getY();
		C[2] = pC->getZ();
		double* const Ap = A;
		double* const Bp = B;
		double* const Cp = C;
		// cout << "trangle points ---------------------" << endl;
		// printArray(Ap,3);
		// printArray(Bp,3);
		// printArray(Cp,3);
		// printArray(ept,3);
		// printArray(dp,3);
		// cout << "------------------------------------" << endl;
		double ax = Ap[0];
		double ay = Ap[1];
		double az = Ap[2];

		double bx = Bp[0];
		double by = Bp[1];
		double bz = Bp[2];

		double cx = Cp[0];
		double cy = Cp[1];
		double cz = Cp[2];

		double lx = ept[0];
		double ly = ept[1];
		double lz = ept[2];

		double dx = dp[0];
		double dy = dp[1];
		double dz = dp[2];
		double M = ((az -cz)*dy-(ay-cy)*dz)*(ax-bx)-((az-cz)*dx-(ax-cx)*dz)*(ay-by) + ((ay-cy)*dx-(ax-cx)*dy)*(az-bz);

		// cout << "z: " << z << endl;
		// calculate Beta if 0 < Beta <= 1 calualte other 
		double M1 = ((az -cz)*dy-(ay-cy)*dz)*(ax-lx)-((az-cz)*dx-(ax-cx)*dz)*(ay-ly)+((ay-cy)*dx-(ax-cx)*dy)*(az-lz);
		double beta = M1/M;
		if (beta >= 0 && beta <= 1){
			
			double M2 = ((az-lz)*dy-(ay-ly)*dz)*(ax-bx)- ((az-lz)*dx-(ax-lx)*dz)*(ay-by)+((ay-ly)*dx - (ax-lx)*dy)*(az-bz);
			// calcualte gama if 0 < gama < 1 && beta + gama <=1
			double gama = M2/M;
			if (gama >=0 && gama <=1 && (gama+beta)<=1){
				// cout << "Beta: " << beta << endl;
				// cout << "gama: " << gama << endl;
				double tM = ((ay-ly)*(az-cz)-(ay-cy)*(az-lz))*(ax-bx)-((ax-lx)*(az-cz)-(ax-cx)*(az-lz))*(ay-by)+((ax-lx)*(ay-cy)-(ax-cx)*(ay-ly))*(az-bz);
				double t = tM/M;
				// cout << "------------------------------------" << endl;
				// cout << "t: " << t <<endl;
				if (t>0){
					if (t < min){
						min = t;
						minFace = run;
					}
					tv.push_back(t);
				}
			}
		} 
	}
	// cout << "t: " << 0 << endl;
	// cout << tv.size()<<endl;
	pair<double, int> p;
	if (min ==  1000){
		p.first = -42;
		p.second = -1;
		return p;
	}
	p.first = min;
	p.second = minFace;
	// cout << "calculate intersection end --------------------------------" << endl;
	return p;
}

pair<double, int> imageCreator::checkSphereIntersection(double* U, double* L){
	// cout << "entering checkSphereIntersection " << endl;
	// printArray(U, 3, "direction");
	// printArray(L, 3, "base");
	vector <pair<double, int>>* p = new vector<pair<double, int>>();
	for (unsigned int i=0; i< spheres->size(); i++){
		sphere* temp = spheres->at(i);
		double* const center = temp->getCenter();
		double radius = temp->getRadius();
		double* const Tv = neo.subtract(center, L);
		double v = neo.dotProduct(Tv,U);
		double bsq = neo.dotProduct(Tv,Tv);
		double disk = (radius*radius) - (bsq - (v*v));
		// cout << "disk " << disk << endl;
		if (disk> 0){
			pair<double, int> match;
			match.first = v-sqrt(disk);
			match.second = i;
			p->push_back(match);
		}

	}
	pair<double, int> minimum;
	minimum.first = 100;
	for (unsigned int j =0; j<p->size(); j++){
		if (p->at(j).first < minimum.first && p->at(j).first > 0.00000001){
			minimum.first = p->at(j).first;
			minimum.second = p->at(j).second;
		}
		// cout << "first: " << minimum.first << endl;
		// cout << "second: " << minimum.second << endl;
	}
	// cout << "exiting checkSphereIntersection" << endl;
	if (minimum.first ==100){
		pair<double, int> empty;
		empty.first = 0;
		empty.second = -1;
		return empty;
	}

	return minimum;

}

void imageCreator::writeOut(){
	// cout << "Target resolution: " << width << " by " << height<< endl;
	// cout << "Polygon Count " << faces->size() << endl;
	// cout << "Depth t runs from " << minD << " to " << maxD << endl;
	outfileName = outfileName.append(".ppm");
	ofstream out(outfileName);
	out <<"P" << 3 <<endl;
	out << width << " " << height << " " << 255 <<endl;
	for (int i =height-1; i>0; i--){
		vector<int*>* temp = newColors->at(i);
		for (int j =0; j<width; j++){
			// cout << "x at: " << j << endl;
			int* currColor = temp->at(j);
			out << currColor[0] << " " << currColor[1] << " " << currColor[2] << " ";
		}
			out << endl;
	}
	out.close();
}

int* imageCreator::ray_trace(double* ray, double* pixpt, double* accum, double* refatt, int level){
	// cout << "entering ray_trace" << endl;
	// printArray(ray, 3, "ray direction");
	// printArray(pixpt, 3, "base point");
	pair<double, int> intersec = checkSphereIntersection(ray, pixpt);
	if (intersec.first > 0.000){
		// cout << "found intersection --------------------------------" << endl;
		double* pt = neo.add(pixpt, neo.scale(intersec.first, ray));
		double * intersectPt = new double[3];
		intersectPt[0] = pt[0];
		intersectPt[1] = pt[1];
		intersectPt[2] = pt[2];
		double* const intersecPt = intersectPt;
		// printArray(intersecPt, 3, "intersection point");

		sphere* s = spheres->at(intersec.second);
		double* snrm = new double[3]; 
		snrm = neo.subtract(intersecPt, s->getCenter());
		// printArray(snrm, 3, "surface Normal");
		// printArray(intersecPt, 3, "intersect");
		// printArray(s->getCenter(), 3, "sphrer center");
		snrm = neo.normalize(snrm);
		// printArray(snrm, 3, "surface Normalized");
		double snrmpt [3];
		snrmpt[0] = snrm[0];
		snrmpt[1] = snrm[1];
		snrmpt[2] = snrm[2];
		double* const bs = snrmpt;

		// printArray(bs, 3, "surface NormalizedPt");
		double* const co = s->getColor();
		double* const am = ambient;
		double* realcolor = new double[3];
		double* tempColor = neo.pair_wise(co, am);
		realcolor[0] = tempColor[0];
		realcolor[1] = tempColor[1];
		realcolor[2] = tempColor[2];
		double* color = realcolor;
		// printArray(color, 3, "starting color");
		// printArray(am, 3, "ambient");
		// printArray(co, 3, "sphere color pt");
		for(unsigned int i =0; i< lights->size(); i++){
			light* lt = lights->at(i);
			double* lpos = lt->getPos();
			double* lcol = lt->getColor();
			// printArray(lpos, 3, "light pos");
			// printArray(lcol, 3, "lcol color");
			// printArray(intersecPt, 3, "Intersect");
			double* toL = new double[3]; 
			toL = neo.subtract(lpos, intersecPt);
			toL = neo.normalize(toL);
			double bss [3];
			bss[0] = toL[0];
			bss[1] = toL[1];
			bss[2] = toL[2];
			double* const toLpt = bss;
			// printArray(toLpt, 3, "toL normalize");
			// cout << "dotProduct: " << neo.dotProduct(bs, toLpt) << endl;
			if (neo.dotProduct(bs, toLpt) > 0){
				// cout << "found intersection----------------------" << endl;
				double scale = neo.dotProduct(bs, toLpt);
				// cout << "scale1: " << scale << endl;
				double* p = neo.pair_wise(kd, lcol);
				// printArray(p, 3, "p array");
				double* f = neo.scale(scale, p);
				// printArray(f, 3, "f array");
				// printArray(color, 3, "color before addion of f---------------");
				color = neo.add(color, f);
				// printArray(color, 3, "color after addion of f ---------------");

				double* toC = neo.subtract(pixpt, intersecPt);
				// printArray(toC, 3, " toC before normalized");
				toC = neo.normalize(toC);
				// printArray(toC, 3, " toC after normalized");

				scale = (2*neo.dotProduct(bs, toLpt));
				// cout << "scale2: " << scale << endl;
				double* spR = neo.scale(scale, bs);
				// printArray(spR, 3, "spR scaled");
				spR = neo.subtract(spR, toLpt);
				// printArray(toLpt,3, "toLpt for reference");
				// printArray(spR, 3, "spR subtacted from toLpt");
				scale = pow(neo.dotProduct(toC, spR), 16);

				// printArray(co, 3, "sphere color pt");
				// printArray(lcol, 3, "light color2");
				double* pw = neo.pair_wise(ks, lcol);
				// printArray(pw, 3, "pw pair of co a lcol");

				double* temp = neo.scale(scale, pw);
				// printArray(temp, 3, "temp scaled of pw");
				// printArray(color, 3, "color befor addtion of temp-----------");
				color = neo.add(color, temp);
				// printArray(color, 3, "color after addtion of temp-----------");
			}
		}
		for (int i =0; i<3; i++){
			accum[i] += refatt[i]*color[i];
		}
		if (level > 0){
			// printArray(ray, 3, "ray--------");
			double* Uinv = neo.scale(-1, ray);
			double* inv = new double[3];
			inv[0] = Uinv[0];
			inv[1] = Uinv[1];
			inv[2] = Uinv[2];
			double* const Uinvpt = inv;
			// printArray(Uinvpt, 3, "inverse ray");
			double Nscale = 2*neo.dotProduct(bs,Uinvpt);
			double* NUinv = neo.scale(Nscale, bs);
			double* refR = neo.subtract(NUinv, Uinvpt); // new ray direction
			refR = neo.normalize(refR);
			// printArray(refR, 3, "new ray direction");
			double* pwRefatt = neo.pair_wise(kr, refatt); // whats refected
			ray_trace(refR, intersecPt, accum, pwRefatt, level-1);
		}

	}
	// cout << "r: " << color[0] << " g: " << color[1] << " b: " << color[2] << "-----------"<< endl;
	int* i = new int[3];
	double r =  max(0.0, 255 * (accum[0])); 
	if (r>255)
		r=255;
	double b = max(0.0, 255 * (accum[2] )); 
	if (b>255)
		b = 255;
	double g = max(0.0, 255 * (accum[1]));
	if (g>255)
		g=255;
	// cout << "r: " << r << " g: " << g << " b: " << b << endl;
	i[0] = (int) r;
	i[1] = (int) g;
	i[2] = (int) b;
	// cout << "exiting ray_trace" << endl;
	return i;
}

int* imageCreator::calculatePolyColor(double* intersecPt, int triPos, double* L){
	// point* f = faces->at(triPos);
	double* snrm = polyHelper(triPos);
	printArray(snrm, 3, "surface Normal");
	printArray(intersecPt, 3, "intersect");
	double snrmpt [3];
	snrmpt[0] = snrm[0];
	snrmpt[1] = snrm[1];
	snrmpt[2] = snrm[2];
	double* const bs = snrmpt;

	// printArray(bs, 3, "surface NormalizedPt");
	double* const co = mat;
	double* const am = ambient;
	double* realcolor = new double[3];
	double* tempColor = neo.pair_wise(co, am);
	realcolor[0] = tempColor[0];
	realcolor[1] = tempColor[1];
	realcolor[2] = tempColor[2];
	double* color = realcolor;
	// printArray(color, 3, "starting color");
	// printArray(am, 3, "ambient");
	// printArray(co, 3, "sphere color pt");
	for(unsigned int i =0; i< lights->size(); i++){
		light* lt = lights->at(i);
		double* lpos = lt->getPos();
		double* lcol = lt->getColor();
		// printArray(lpos, 3, "light pos");
		// printArray(lcol, 3, "lcol color");
		// printArray(intersecPt, 3, "Intersect");
		double* toL = new double[3]; 
		toL = neo.subtract(lpos, intersecPt);
		toL = neo.normalize(toL);
		double bss [3];
		bss[0] = toL[0];
		bss[1] = toL[1];
		bss[2] = toL[2];
		double* const toLpt = bss;
		// printArray(toLpt, 3, "toL normalize");
		// cout << "dotProduct: " << neo.dotProduct(bs, toLpt) << endl;
		if (neo.dotProduct(bs, toLpt) > 0){
			// cout << "found intersection----------------------" << endl;
			double scale = neo.dotProduct(bs, toLpt);
			// cout << "scale1: " << scale << endl;
			double* p = neo.pair_wise(kd, lcol);
			// printArray(p, 3, "p array");
			double* f = neo.scale(scale, p);
			// printArray(f, 3, "f array");
			// printArray(color, 3, "color before addion of f---------------");
			color = neo.add(color, f);
			// printArray(color, 3, "color after addion of f ---------------");

			double* toC = neo.subtract(L, intersecPt);
			// printArray(toC, 3, " toC before normalized");
			toC = neo.normalize(toC);
			// printArray(toC, 3, " toC after normalized");

			scale = (2*neo.dotProduct(bs, toLpt));
			// cout << "scale2: " << scale << endl;
			double* spR = neo.scale(scale, bs);
			// printArray(spR, 3, "spR scaled");
			spR = neo.subtract(spR, toLpt);
			// printArray(toLpt,3, "toLpt for reference");
			// printArray(spR, 3, "spR subtacted from toLpt");

			// cout << "scale2: " << scale << endl;
			// cout << " toc spr dotProduct: " << neo.dotProduct(toC, spR) << endl;
			scale = pow(neo.dotProduct(toC, spR), 16);
			// cout << "scale3: " << scale << endl;

			// printArray(co, 3, "sphere color pt");
			// printArray(lcol, 3, "light color2");
			double* pw = neo.pair_wise(ks, lcol);
			// printArray(pw, 3, "pw pair of co a lcol");

			double* temp = neo.scale(scale, pw);
			// printArray(temp, 3, "temp scaled of pw");
			// printArray(color, 3, "color befor addtion of temp-----------");
			color = neo.add(color, temp);
			// printArray(color, 3, "color after addtion of temp-----------");
		}
	}
	// cout << "r: " << color[0] << " g: " << color[1] << " b: " << color[2] << "-----------"<< endl;
	int* i = new int[3];
	double r =  max(0.0, 255 * (color[0])); 
	if (r>255)
		r=255;
	double b = max(0.0, 255 * (color[2])); 
	if (b>255)
		b = 255;
	double g = max(0.0, 255 * (color[1]));
	if (g>255)
		g=255;
	// cout << "r: " << r << " g: " << g << " b: " << b << endl;
	i[0] = (int) r;
	i[1] = (int) g;
	i[2] = (int) b;
	return i;
}

double* imageCreator::polyHelper(int triPos){
	point* const triangle = faces->at(triPos);
		double pointA = triangle->getX();
		double pointB = triangle->getY();
		double pointC = triangle->getZ();

		point* const pA = points->at(pointA);
		double A [3];
		A[0] = pA->getX();
		A[1] = pA->getY();
		A[2] = pA->getZ();
		point* const pB = points->at(pointB);
		double B [3];
		B[0] = pB->getX();
		B[1] = pB->getY();
		B[2] = pB->getZ();
		point* const pC = points->at(pointC);
		double C[3];
		C[0] = pC->getX();
		C[1] = pC->getY();
		C[2] = pC->getZ();
		double * aToB = neo.subtract(A, B);
		aToB = neo.normalize(aToB);
		double* aToC = neo.subtract(A, C);
		aToC = neo.normalize(aToC);
		double* snrm = neo.crossProduct(aToB, aToC);
		return snrm;
}

int* imageCreator::calculateColor(double t){
	static int color [3];
	double ratio = 2 * (t - minD) / (maxD - minD);
	// cout << "ratio: " << ratio << endl;
	int r = max(0.0, 255 * (1 - ratio)); 
	int b = max(0.0, 255 * (ratio - 1)); 
	int g = (b - r);
	if (g > 255)
		g = 255;
	color[0] = r;
	color[1] = g;
	color[2] = b;
	// cout << "r: " << r << " g: " << g << " b: " << b << endl;
	return color;
}

void imageCreator::printCameraCoord(){
	cout << "Camera Coordinate System --------------" << endl;
	cout << "W " << WC[0] << " ," << WC[1] << " ," << WC[2] << endl;
	cout << "U " << UC[0] << " ," << UC[1] << " ," << UC[2] << endl;
	cout << "V " << VC[0] << " ," << VC[1] << " ," << VC[2] << endl;
}

void imageCreator::printArray(double* p, int size, string s){
	cout << "printing array " << s <<endl;
	for (int i=0; i<size; i++){
		cout << p[i] << " ";
	}
	cout <<endl;
}

void imageCreator::freeCamerCoord(){
	free(WC);
	free(UC);
	free(VC);
}



void imageCreator::printVector(vector<double>* v){
	cout << "printing vector ---------------------------------- " << endl;
	for (unsigned int i =0; i< v->size(); i++){
		cout << v->at(i) << " ";
	}
	cout << endl;
}