#pragma once
class Matrix;

class Vector {
public:
	std::string Name;
	std::vector<double> Data;

	Vector();
	Vector(int cnt,...);
	Vector(std::vector<double>& _data);

	Vector operator+(Vector& v);
	Vector operator-(Vector& v);
	Vector& operator+=(Vector& v);
	Vector& operator-=(Vector& v);
	
	Vector operator*(double scale);
	Vector operator/(double scale);
	Vector& operator*=(double scale);
	Vector& operator/=(double scale);

	bool isZero();
	bool operator ==(Vector& v);

	double Norm();
	Vector Normal();
	Vector& Normalized();

	static double Dot(Vector& v1, Vector& v2);
	static Vector Cross(Vector& v1, Vector& v2);
	static double Com(Vector& v1, Vector& v2);
	static Vector Proj(Vector& v1, Vector& v2);
	static double TriArea(Vector& v1, Vector& v2);
	static bool isParallel(Vector& v1, Vector& v2);
	static bool isOrthognal(Vector& v1, Vector& v2);
	static double Radians(Vector& v1, Vector& v2);
	static double Angle(Vector& v1, Vector& v2);
	static Vector PlaneNormal(Vector& v1, Vector& v2);
	static bool IsLI(Vector& v1, Vector& v2);

	static std::vector<Vector> Ob(std::vector<Vector>& vs);

};
