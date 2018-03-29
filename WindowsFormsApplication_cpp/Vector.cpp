#include "Matrix.h"
#include "Vector.h"

#include<Windows.h>
#include <iostream>
#include <math.h>

using namespace std;

Vector::Vector() {

}

Vector::Vector(vector<double>& _data) {
	this-> Data = _data;
}

Vector::Vector(int cnt, ...) {
	va_list args;
	va_start(args, cnt);
	for (int i = 0; i < cnt; i++) {
		this->Data.push_back(va_arg(args, double));
	}
	va_end(args);
}

bool Vector::isZero() {
	for (auto& d : this->Data) {
		if (d != 0)
			return false;
	}
	return true;
}

Vector& Vector::operator+=(Vector& v) {
	int len = this->Data.size() > v.Data.size() ? v.Data.size() : this->Data.size();
	for (int i = 0; i < len; i++) {
		this->Data[i] += v.Data[i];
	}
	return *this;
}
Vector Vector::operator+(Vector& v) {
	Vector result = *this;
	return result += v;
}
Vector& Vector::operator-=(Vector& v) {
	int len = this->Data.size() > v.Data.size() ? v.Data.size() : this->Data.size();
	for (int i = 0; i < len; i++) {
		this->Data[i] -= v.Data[i];
	}
	return *this;
}
Vector Vector::operator-(Vector& v) {
	Vector result = *this;
	return result -= v;
}

bool Vector::operator ==(Vector& v) {
	if (this->Data.size() != v.Data.size())
		return false;
	for (int i = 0; i < this->Data.size(); i++) {
		if (this->Data[i] != v.Data[i]) {
			return false;
		}
	}
	return true;
}

Vector& Vector::operator*=(double scale) {
	for (int i = 0; i < this->Data.size(); i++) {
		this->Data[i] *= scale;
	}
	return *this;
}
Vector Vector::operator*(double scale) {
	Vector result = *this;
	return result *= scale;
}
Vector& Vector::operator/=(double scale) {
	for (int i = 0; i < this->Data.size(); i++) {
		this->Data[i] /= scale;
	}
	return *this;
}

Vector Vector::operator/(double scale) {
	Vector result = *this;
	return result /= scale;
}

double Vector::Norm() {
	double result = 0;
	for (auto& d : this->Data) {
		result += d*d;
	}
	result = sqrt(result);
	return result;
}
Vector Vector::Normal() {
	Vector result = *this;
	result /= this->Norm();
	return result;
}
Vector& Vector::Normalized() {
	*this /= this->Norm();
	return *this;
}

double Vector::Dot(Vector& v1, Vector& v2) {
	double result = 0;
	for (int i = 0; i < v1.Data.size() && i < v2.Data.size(); i++) {
		result += v1.Data[i] * v2.Data[i];
	}
	return result;
}

// only for 3d
Vector Vector::Cross(Vector& v1, Vector& v2) {
	Vector result;
	result.Data.push_back(v1.Data[1]*v2.Data[2] - v1.Data[2] * v2.Data[1]);
	result.Data.push_back(-v1.Data[0] * v2.Data[2] + v1.Data[2] * v2.Data[0]);
	result.Data.push_back(v1.Data[0] * v2.Data[1] - v1.Data[1] * v2.Data[0]);
	return result;
}
double Vector::Com(Vector& v1, Vector& v2) {
	double result;
	result = Dot(v1, v2);
	result /= v2.Norm();
	return result;
}
Vector Vector::Proj(Vector& v1, Vector& v2) {
	Vector result = v2;
	result *= Com(v1, v2);
	result /= v2.Norm();
	return result;
}
double Vector::TriArea(Vector& v1, Vector& v2) {
	double result;
	result = v1.Norm()*v2.Norm() / 2.0;
	return result;
}
bool Vector::isParallel(Vector& v1, Vector& v2){
	double r = Radians(v1, v2);
	return (r == 0.0 || r == 1.0 || r == -1.0);
}
bool Vector::isOrthognal(Vector& v1, Vector& v2) {
	double r = Radians(v1, v2);
	return (r == 0.5 || r == -0.5);
}
double Vector::Radians(Vector& v1, Vector& v2) {
	double result;
	result = Dot(v1, v2);
	result /= v1.Norm()*v2.Norm();
	result = acos(result);
	return 0;
}

double Vector::Angle(Vector& v1, Vector& v2) {
	double result;
	result = Radians(v1, v2);
	result *= 180.0 / 3.14159265;
	return 0;
}

Vector Vector::PlaneNormal(Vector& v1, Vector& v2) {
	return Cross(v1,v2).Normalized();
}
bool Vector::IsLI(Vector& v1, Vector& v2) {
	bool result;
	result = v1.Data[0] * v2.Data[1] - v1.Data[1] * v2.Data[0];
	return result;
}

vector<Vector> Vector::Ob(vector<Vector>& vs) {
	vector<Vector> result;
	vector<double> tmpDots;
	for (int i = 0; i < vs.size(); i++) {
		result.push_back(vs[i]);
		for (int j = 0; j < i; j++) {
			result[i] -= result[j]*(Dot(result[i], result[j])/tmpDots[j]);
		}
		tmpDots.push_back(Dot(result[i], result[i]));
	}
	for (auto& v : result) {
		v.Normalized();
	}
	return result;
}