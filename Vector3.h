#pragma once

class vector3{
public:
	double x, y, z;
	vector3(double in_x = 0, double in_y = 0, double in_z = 0){
		x = in_x; y = in_y; z = in_z;
	}
	~vector3(){}

	vector3 operator + (const vector3 &v) const { return vector3(x + v.x, y + v.y, z + v.z); }
	vector3 operator - (const vector3 &v) const { return vector3(x - v.x, y - v.y, z - v.z); }
	vector3 operator + (double v) const { return vector3(x + v, y + v, z + v); }
	vector3 operator - (double v) const { return vector3(x - v, y - v, z - v); }
	vector3 operator * (double v) const { return vector3(x * v, y * v, z * v); }
	vector3 operator % (vector3 &v) { return vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	vector3 normllize() { return  vector3(x / sqrt(x * x + y * y + z * z), y / sqrt(x * x + y * y + z * z), z / sqrt(x * x + y * y + z * z)); }
	vector3 mul(const vector3 &v) const { return vector3(x * v.x, y * v.y, z * v.z); }
	double dot(const vector3 &v) const { return x * v.x + y * v.y + z * v.z; }
};
