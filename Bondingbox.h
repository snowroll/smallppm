#pragma once
#include "Vector3.h"

#define MAX(x, y) ((x > y) ? x : y)

class AABB{
public:
	AABB(){}
	~AABB(){}

	vector3 min, max;  //OBB求交
	void fit(const vector3 &pos){
		if(pos.x < min.x) min.x = pos.x;
		if(pos.y < min.y) min.y = pos.y;
		if(pos.z < min.z) min.z = pos.z;
		max.x = MAX(pos.x, max.x);
		max.y = MAX(pos.y, max.y);
		max.z = MAX(pos.z, max.z);
	}
	void reset(){
		min = vector3(1e20, 1e20, 1e20);
		max = vector3(-1e20, -1e20, -1e20);
	} 
};