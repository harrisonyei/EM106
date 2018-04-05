#include"Matrix.h"
#include"Vector.h"
#include <iostream>
using namespace std;

#ifdef DEBUG_M

bool FloatEqual(double a,double b,double precision = 1e-5){
	return (abs(a-b)<precision);
}

void Matrix::PrintM(Matrix& m) {
	cout << "Name : " << endl;
	for (auto& v : m.Data) {
		for (auto& d : v) {
			cout << d << "\t";
		}
		cout << endl;
	}
	cout << endl;
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

	for (int i = 0; i < m1RowSize; i++) {
		vector<double> row;
		for (int j = 0; j < m2ColSize; j++) {
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

	for (int i = 0; i < m1RowSize; i++) {
		vector<double> row;
		for (int j = 0; j < m2ColSize; j++) {
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
Matrix operator*(Vector& v,Matrix& m){
	Matrix result;
	if(m.Data.size() != v.Data.size())
		throw "invalid";

	int m1RowSize = 1;
	int m2RowSize = m.Data.size();

	int m1ColSize = v.Data.size();
	int m2ColSize = m.Data.front().size();

	for(int i = 0; i < m1RowSize; i++){
		vector<double> row;
		for(int j = 0; j < m2ColSize; j++){
			double value = 0;
			for(int k = 0; k < m1ColSize; k++){
				value += m.Data[k][i] * v.Data[k];
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
double PolyFunction(double* func,int c,double x){
	double y = 0;
	for(int i = 0; i < c; i++){
		double tmp = func[i];
		for(int j = i+1;j < c;j++){
			tmp *= x;
		}
		y += tmp;
	}
	return y;
}
vector<double> NewtonRoot(double* func,double* dfunc,double f,double precision = 1e-5,int maxtry = 10000){
	vector<double> result;
	vector<double> x_start;

	double delta = dfunc[1] * dfunc[1] - 4 * dfunc[0] * dfunc[2];
	double sqrtDelta = sqrt(delta);
	
	double x1 = ((-dfunc[1] + sqrtDelta) / (2 * dfunc[0]));
	double x2 = ((-dfunc[1] - sqrtDelta) / (2 * dfunc[0]));
	
	double y1 = PolyFunction(func,4,x1);
	double y2 = PolyFunction(func,4,x2);

	if(FloatEqual(f,0)){
		bool y10 = FloatEqual(y1,0);
		bool y20 = FloatEqual(y2,0);
		if(y10&&y20){
			result.push_back(x1);
		} else if(y10){
			result.push_back(x1);
			x_start.push_back(x2 + (x2 - x1));
		} else{
			result.push_back(x2);
			x_start.push_back(x1 + (x1 - x2));
		}
	} else if(f > 0){
		if(y1 > y2){
			x_start.push_back(x2 + (x2 - x1));
		} else{
			x_start.push_back(x1 + (x1 - x2));
		}
	} else if(f < 0){
		x_start.push_back(x2 + (x2 - x1));
		x_start.push_back(x1 + (x1 - x2));
		x_start.push_back((x1 + x2) / 2);
	}

	double sol,m,b;

	for(auto& x_sol : x_start){
		double tmp_x_sol = 0;
		sol = PolyFunction(func,4,x_sol);
		while(abs(sol)>precision || abs(x_sol - tmp_x_sol)>precision){
			tmp_x_sol = x_sol;
			m = PolyFunction(dfunc,3,x_sol);
			b = sol - m * x_sol;
			x_sol = b / -m;
			sol = PolyFunction(func,4,x_sol);
			if(!(maxtry--)){
				throw "Fail to find answer";
			}
		}
		result.push_back(x_sol);
	}

	return result;
}

void Matrix::Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue) {
	int rank = m.Data.size();
	vector<double> lamda;
	// rank = 3
	if(rank == 3){
		double d = (m.Data[0][0] * m.Data[1][1] * m.Data[2][2] + m.Data[1][0] * m.Data[2][1] * m.Data[0][2] + m.Data[0][1] * m.Data[1][2] * m.Data[2][0] \
			- m.Data[0][2] * m.Data[1][1] * m.Data[2][0] - m.Data[1][2] * m.Data[2][1] * m.Data[0][0] - m.Data[0][1] * m.Data[1][0] * m.Data[2][2]);
		double c = (-m.Data[0][0] * m.Data[1][1] - m.Data[1][1] * m.Data[2][2] - m.Data[2][2] * m.Data[0][0]\
			+ m.Data[0][1] * m.Data[1][0] + m.Data[1][2] * m.Data[2][1] + m.Data[0][2] * m.Data[2][0]);
		double b = (m.Data[0][0] + m.Data[1][1] + m.Data[2][2]);
		double a = (-1);
		
		double t1 = (36 * a*b*c - 8 * b*b*b - 108 * a*a*d);
		double t2 = (12 * a*c - 4 * b*b);
		double f = t1*t1 + t2*t2*t2;

		double func[] = {a,b,c,d};
		double dfunc[] = {3*a,2*b,c};

		lamda = (NewtonRoot(func,dfunc,f,1e-8));
	} else if(rank == 2){
		double c = (m.Data[0][0] * m.Data[1][1] - m.Data[0][1] * m.Data[1][0]);
		double b = (-m.Data[0][0] - m.Data[1][1]);
		double a = (1);

		double f = b*b - 4 * a*c;
		if(f == 0){
			lamda.push_back(-b/(2*a));
		} else if(f > 0){
			double sqrtf = sqrt(f);
			lamda.push_back((-b + sqrtf) / (2 * a));
			lamda.push_back((-b - sqrtf) / (2 * a));
		} else{
			throw "no ans";
		}
	} else{
		double b = (m.Data[0][0]);
		double a = (-1);
		lamda.push_back(b);
	}
	for(auto& l : lamda){
		Matrix tmpM = m;
		for(int i = 0;i < tmpM.Data.size();i++){
			tmpM.Data[i][i] -= l;
		}
		cout << l << endl;
		G_Eliminate(tmpM,true);
		PrintM(tmpM);
	}
	// not finished
}
void Matrix::PM_Eigen(Matrix& m, Matrix& eigenVector, Matrix& eigenValue) {
	// not finished
	double current_lamda=1,pre_lamda=0;
	Vector u;
	Matrix Au;
	u.Data = {1,0,0};
	while(abs(current_lamda-pre_lamda)>1e-5){
		Au = m*u;
		u.Data = Au.T().Data[0];
		u.Normalized();
		pre_lamda = current_lamda;
		current_lamda = (u*Au).Data[0][0];
	}
	vector<vector<double>> vs;
	vs.push_back(u.Data);
	vs.push_back((m*u).T().Data[0]);
	vs.push_back((m*Vector(vs.back())).T().Data[0]);
	for(auto& v : vs){
		for(auto& i : v){
			cout << i << " ";
		}
		cout << endl;
	}
	cout << current_lamda << endl;
}

void CancelAt(vector<double>& src ,vector<double>& tar,int idx) {
	if (src[idx] != 0 && tar[idx] != 0) {
		double diff;
		if(FloatEqual(tar[idx],src[idx])){
			diff = 1;
		} else{
			diff = tar[idx] / src[idx];
		}
		tar[idx] = 0;
		for (int i = idx+1; i < src.size(); i++) {
			tar[i] -= diff*src[i];
			if(FloatEqual(tar[i],0)){
				tar[i] = 0;
			}
		}
	}
}

// need debug for {1,0,0}
//                {0,1,0}
//                {0,0,0}
std::vector<double> Matrix::SolveLinearSys(Matrix& m1, Matrix& m2) {
	// not finished
	vector<double> result;
	for (int i = 0; i < m1.Data.size(); i++) {
		m1.Data[i].push_back(m2.Data[i][0]);
	}
	G_Eliminate(m1,true);
	for (int i = 0; i < m1.Data.size(); i++) {
		if(m1.Data[i][i] == 0){
			if(m1.Data[i].back() == 0)
				throw "inf Answers";
			else
				throw "No Answer";
		}
		result.push_back(m1.Data[i].back());
	}
	return result;
}

// need debug for {1,0,0}
//                {0,1,0}
//                {0,0,0}
Vector Matrix::LeastSquare(Matrix& m, Matrix& v) {
	Matrix m1,m2;
	for(int i = 0; i < m.Data.front().size();i++){
		vector<double> row;
		for(int j = 0; j < m.Data.front().size();j++){
			row.push_back(0);
			for(int k = 0; k < m.Data.size();k++){
				row.back() += 2 * m.Data[k][j] * m.Data[k][i];
			}
		}
		m1.Data.push_back(row);
		m2.Data.push_back({0});
		for(int k = 0; k < v.Data.size();k++){
			m2.Data.back()[0] += 2 * v.Data[k][0] * m.Data[k][i];
		}
	}
	return Vector(SolveLinearSys(m1,m2));
}

std::vector<Matrix> Matrix::rref(Matrix& m) {
	vector<Matrix> result;
	result.push_back(G_Eliminate(m,true));
	result.push_back(G_Eliminate(m,false));
	return result;
}

Matrix& Matrix::G_Eliminate(Matrix& m,bool up){
	int idx = 0;
	for(int i = 0; i < m.Data.size(); i++){
		for(int end = m.Data.size() - 1; !m.Data[i][idx] && end > i; end--){
			vector<double> tmp;
			tmp = m.Data[i];
			m.Data[i] = m.Data[end];
			m.Data[end] = tmp;
		}

		for(int j = i + 1; j < m.Data.size(); j++){
			CancelAt(m.Data[i],m.Data[j],idx);
		}

		if(m.Data[i][idx] != 0){
			for(int k = idx + 1; k < m.Data[i].size(); k++){
				m.Data[i][k] /= m.Data[i][idx];
			}
			m.Data[i][idx] = 1;
		}
		idx += 1;
	}
	for(int i = m.Data.size() - 1; i >= 0; i--){
		idx -= 1;
		for(int j = i - 1; j >= 0; j--){
			CancelAt(m.Data[i],m.Data[j],idx);
		}
	}
	if(!up){
		for(int i = 0,j = m.Data.size() - 1; i < j;i++,j--){
			vector<double> tmp;
			tmp = m.Data[i];
			m.Data[i] = m.Data[j];
			m.Data[j] = tmp;
		}
	}
	return m;
}
