#pragma once
#include<string>
#include<vector>
#define DEBUG_M

class Vector;

class Matrix {
public:
	std::string Name;
	std::vector<std::vector<double>> Data;

	Matrix operator+(Matrix& v);
	Matrix operator-(Matrix& v);
	Matrix operator*(Matrix& m);
	Matrix operator*(Vector& v);
	Matrix operator*(double scale);
	Matrix operator/(double scale);

	friend Matrix operator*(Vector& v,Matrix& m);

	Matrix& operator+=(Matrix& v);
	Matrix& operator-=(Matrix& v);
	Matrix& operator*=(Matrix& m);
	Matrix& operator*=(Vector& v);
	Matrix& operator*=(double scale);
	Matrix& operator/=(double scale);

	Matrix T();
	double Det();
	Matrix Inv();
	Matrix Adj();

	static Matrix& Trans(Matrix& m);
	static int Rank(Matrix& m);
	static Matrix& Inverse(Matrix& m);
	static Matrix& Adjoint(Matrix& m);
	static void Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue);
	static void PM_Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue);
	static std::vector<double> SolveLinearSys(Matrix& m1, Matrix& m2);
	static Vector LeastSquare(Matrix& m, Matrix& v);
	static std::vector<Matrix> rref(Matrix& m);
	static Matrix& G_Eliminate(Matrix& m,bool up);
#ifdef DEBUG_M
	static void PrintM(Matrix& m);
#endif // DEBUG_M

private:
	Matrix ShrinkMatrix(int row, int col);
};