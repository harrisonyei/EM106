#include"Matrix.h"
#include"Vector.h"
#include <iostream>
using namespace std;

#ifdef DEBUG_M
void Matrix::PrintM(Matrix& m) {
	for (auto& v : m.Data) {
		for (auto& d : v) {
			cout << d << " ";
		}
		cout << endl;
	}
}
#endif // DEBUG_M

Matrix Matrix::operator+(Matrix& v) {
	Matrix result = *this;
	result += v;
	return result;
}
Matrix Matrix::operator-(Matrix& v) {
	Matrix result = *this;
	result -= v;
	return result;
}
Matrix Matrix::operator*(Matrix& m) {
	Matrix result;
	if (this->Data.begin()->size() != m.Data.size())
		throw "invalid";
	int m1RowSize = this->Data.size();
	int m2RowSize = m.Data.size();

	int m1ColSize = this->Data.begin()->size();
	int m2ColSize = m.Data.begin()->size();

	for (int i = 0; i < m2RowSize; i++) {
		vector<double> row;
		for (int j = 0; j < m1ColSize; j++) {
			double value = 0;
			for (int k = 0; k < m1ColSize; k++) {
				value += this->Data[i][k] * m.Data[k][j];
			}
			row.push_back(value);
		}
		result.Data.push_back(row);
	}
	return result;
}
Matrix Matrix::operator*(Vector& v) {
	Matrix result;
	if (this->Data.begin()->size() != v.Data.size())
		throw "invalid";

	int m1RowSize = this->Data.size();
	int m2RowSize = v.Data.size();

	int m1ColSize = this->Data.begin()->size();
	int m2ColSize = 1;

	for (int i = 0; i < m2RowSize; i++) {
		vector<double> row;
		for (int j = 0; j < m1ColSize; j++) {
			double value = 0;
			for (int k = 0; k < m1ColSize; k++) {
				value += this->Data[i][k] * v.Data[k];
			}
			row.push_back(value);
		}
		result.Data.push_back(row);
	}
	return result;
}
Matrix Matrix::operator*(double scale) {
	Matrix result = *this;
	result *= scale;
	return result;
}
Matrix Matrix::operator/(double scale) {
	Matrix result = *this;
	result /= scale;
	return result;
}
Matrix& Matrix::operator+=(Matrix& v) {
	if (this->Data.size() != v.Data.size() || this->Data.begin()->size() != v.Data.begin()->size())
		throw "invalid";
	for (int i = 0; i < this->Data.size(); i++) {
		for (int j = 0; j < this->Data.begin()->size(); j++) {
			this->Data[i][j] += v.Data[i][j];
		}
	}
	return *this;
}
Matrix& Matrix::operator-=(Matrix& v) {
	if (this->Data.size() != v.Data.size() || this->Data.begin()->size() != v.Data.begin()->size())
		throw "invalid";
	for (int i = 0; i < this->Data.size(); i++) {
		for (int j = 0; j < this->Data.begin()->size(); j++) {
			this->Data[i][j] -= v.Data[i][j];
		}
	}
	return *this;
}
Matrix& Matrix::operator*=(Matrix& m) {
	*this = *this*m;
	return *this;
}
Matrix& Matrix::operator*=(Vector& v) {
	*this = *this*v;
	return *this;
}
Matrix& Matrix::operator*=(double scale) {
	for (int i = 0; i < this->Data.size(); i++) {
		for (int j = 0; j < this->Data.begin()->size(); j++) {
			this->Data[i][j] *= scale;
		}
	}
	return *this;
}
Matrix& Matrix::operator/=(double scale) {
	for (int i = 0; i < this->Data.size(); i++) {
		for (int j = 0; j < this->Data.begin()->size(); j++) {
			this->Data[i][j] /= scale;
		}
	}
	return *this;
}

Matrix Matrix::T() {
	Matrix result;
	for (int j = 0; j < this->Data.begin()->size(); j++) {
		vector<double> row;
		for (int i = 0; i < this->Data.size(); i++) {
			row.push_back(this->Data[i][j]);
		}
		result.Data.push_back(row);
	}
	return result;
}
double Matrix::Det() {
	double result = 0;
	if (this->Data.size() == 1) {
		result = this->Data[0][0];
	}
	else if (this->Data.size() == 2) {
		result = this->Data[0][0] * this->Data[1][1] - this->Data[0][1] * this->Data[1][0];
	}
	else if (this->Data.size() == 3) {
		result = this->Data[0][0] * this->Data[1][1]* this->Data[2][2] + this->Data[1][0] * this->Data[2][1] * this->Data[0][2] + this->Data[0][1] * this->Data[1][2] * this->Data[2][0] \
			- this->Data[0][2] * this->Data[1][1] * this->Data[2][0] - this->Data[1][2] * this->Data[2][1] * this->Data[0][0] - this->Data[0][1] * this->Data[1][0] * this->Data[2][2];
	}
	else {
		for (int i = 0; i < this->Data.begin()->size(); i++) {
			result += ((i & 1) ? -this->Data[0][i] : this->Data[0][i])*ShrinkMatrix(0, i).Det();
		}
	}
	return result;
}

Matrix Matrix::ShrinkMatrix(int row, int col) {
	Matrix result;
	for (int i = 0; i < this->Data.size(); i++) {
		if (i != row) {
			vector<double> row;
			for (int j = 0; j < this->Data.begin()->size(); j++) {
				if (j != col) {
					row.push_back(this->Data[i][j]);
				}
			}
			result.Data.push_back(row);
		}
	}
	return result;
}

Matrix Matrix::Inv() {
	double det = this->Det();
	Matrix result;

	if (det == 0) {
		throw "Error";
	}

	result = this->Adj();
	result /= det;

	return result;
}
Matrix Matrix::Adj() {
	Matrix result;
	for (int j = 0; j < this->Data.size(); j++) {
		vector<double> row;
		for (int i = 0; i < this->Data.size(); i++) {
			row.push_back((((i & 1) ^ (j & 1)) ? -1 : 1)*ShrinkMatrix(i, j).Det());
		}
		result.Data.push_back(row);
	}
	return result;
}

Matrix& Matrix::Trans(Matrix& m) {
	m = m.T();
	return m;
}
int Matrix::Rank(Matrix& m) {
	int result = 0;
	vector<Vector> vs;
	for (auto& v : m.Data) {
		vs.push_back(v);
	}
	vs = Vector::Ob(vs);
	for (auto& v : vs) {
		if (!v.isZero()) {
			result += 1;
		}
	}
	return result;
}
Matrix& Matrix::Inverse(Matrix& m) {
	m = m.Inv();
	return m;
}
Matrix& Matrix::Adjoint(Matrix& m) {
	m = m.Adj();
	return m;
}
void Matrix::Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue) {
	// not finished
}
void Matrix::PM_Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue) {
	// not finished
}

void CancelAt(vector<double>& src ,vector<double>& tar,int idx) {
	if (src[idx] != 0 && tar[idx] != 0) {
		double diff = tar[idx] / src[idx];
		tar[idx] = 0;
		for (int i = idx+1; i < src.size(); i++) {
			tar[i] -= diff*src[i];
		}
	}
}

std::vector<double> Matrix::SolveLinearSys(Matrix& m1, Matrix& m2) {
	// not finished
	vector<double> result;
	int idx = 0;
	for (int i = 0; i < m1.Data.size(); i++) {
		m1.Data[i].push_back(m2.Data[i][0]);
	}
	for (int i = 0; i < m1.Data.size(); i++) {

		for (int end = m1.Data.size() - 1; !m1.Data[i][idx] && end > i; end--) {
			vector<double> tmp;
			tmp = m1.Data[i];
			m1.Data[i] = m1.Data[end];
			m1.Data[end] = tmp;
		}

		for (int j = i+1; j < m1.Data.size(); j++) {
			CancelAt(m1.Data[i], m1.Data[j], idx);
		}
		idx += 1;
	}
	for (int i = m1.Data.size()-1; i >= 0; i--) {
		idx -= 1;
		for (int j = i - 1; j >= 0; j--) {
			CancelAt(m1.Data[i], m1.Data[j], idx);
		}
	}
	PrintM(m1);
	for (int i = 0; i < m1.Data.size(); i++) {
		result.push_back(m1.Data[i].back()/m1.Data[i][i]);
	}
	return result;
}
Vector Matrix::LeastSquare(Matrix& m, Matrix& v) {
	// not finished
	return Vector();
}
std::vector<Matrix> Matrix::rref(Matrix& m) {
	// not finished
	return std::vector<Matrix>();
}