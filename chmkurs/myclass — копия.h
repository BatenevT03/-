
#define M 8
//----------------------------
//#include "local_matrix.h"
//#include <alloc.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ranges>
#include <iomanip>
#include <chrono>  
using namespace std;
//----------------------------
double right(vector<double>&x);
void mul_c_vector(vector<double>&a, vector<double>&b);

// Локальные области
struct local //i-ый КЭ
{
	static int count;
	vector<int> mas;
	//int mas[8];//номера узлов в КЭ
	double lambda;// лямда для КЭ
};

//Координаты узлов
struct cross
{
	vector<double> mas;
	//double mas[3];
	static int count;
};



class MATRIX
{
public:
	int n;
	vector<double>d;
	vector<double>gg;
	vector<int>ig;
	vector<int>jg;


	void solution_x_l(vector<double>& f, vector<double>& x);
	void solution_x_u(vector<double>& f, vector<double>& x);
	void solution_x(vector<double>& f, vector<double>& x);
	void mul_matrix_vector(vector<double>& vec, vector<double>& result);
	void SST(MATRIX *);
	~MATRIX();
};



class GLOBAL_MATRIX : public MATRIX
{

	


public:
	local* matr;//i-ый КЭ
	cross* set;
	vector<double>f;
	vector<double>x;
	vector<double>vect_v;
	GLOBAL_MATRIX();
	~GLOBAL_MATRIX();
	void UpdateGaussPointValues( vector<double>& solution);
	void read_local();
	void read_cross();
	void formier_profil();
	void formier_matrix();
	void add(int,int,double);
	void MSG();
	void kraev_1();
	void kraev_2();
	void generate();
	double U(double nx, double y, double z);
	double get_psi(int num_fe, int num_basis, double x,double y, double z);
	const char *FILE_LOCAL_NUMBER;
	const char *FILE_CROSS;
	const char *FILE_OUT;
	const char *FILE_CRAEV_1;
	const char *FILE_CRAEV_2;
	double POR;
	void print_result();
};

class Interpolation : public GLOBAL_MATRIX
{
private:
	vector<double> x, y, z,f,phi,psi;
public:
	Interpolation()
	{
		double a, b, c, d;
		ifstream in("out.txt");
		while (!in.eof())
		{
			in >> a >> b >> c >> d;
			x.push_back(a);
			y.push_back(b);
			z.push_back(c);
			f.push_back(d);
		}
	}
	double U(double xi, double eta, double theta);
};


struct PointData {
	double x, y, z;  // Координаты точки
	double value;     // Значение решения в точке
};

