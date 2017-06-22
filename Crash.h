#pragma once 
#include "Vector3.h"

class Crash_point{
public:
	Crash_point(){}
	~Crash_point(){}

	vector3 f, pos, nrm, flux;
	double r2;
	unsigned int n;  
	int pix;
};

class List{
public:
	List(){}
	~List(){}

	Crash_point *id;
	List *next;
};
