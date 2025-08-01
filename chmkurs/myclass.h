
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
	//local* matr;//i-ый КЭ
	//cross* set;
	vector<local> matr;
	vector<cross> set;
	vector<double>f;
	vector<double>x;
	//vector<double>vect_v;
	GLOBAL_MATRIX();
	
	~GLOBAL_MATRIX();
	void read_local(string);
	void read_cross(string);
	void formier_profil();
	void formier_matrix();
	void add(int,int,double);
	void MSG();
	void kraev_1();
	void kraev_2();
	void generate();
	int find6el(double x, double y, double z);
	
	void Method_Newton(int nf1, vector<double>& a,bool sign);

	void SetF();
	double U6(double x, double y, double z);
	//void Ugauss(int a);
	double U(double nx, double y, double z);
	double get_psi(int num_fe, int num_basis, double x,double y, double z);
	const char *FILE_LOCAL_NUMBER;
	const char *FILE_CROSS;
	const char *FILE_OUT;
	const char *FILE_CRAEV_1;
	const char *FILE_CRAEV_2;
	double POR;
	void print_result();
	void Receiver(ofstream& out, double len, double start_x, double start_z);
	void print_receiver();
	//void Interpolation();
	GLOBAL_MATRIX(vector<cross>& Cross, vector<local>& Matr, vector<double>& q);
};




struct PointData
{
	double x;
	double y;
	double z;
	double value;
	PointData(double x, double y, double z, double value)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->value = value;
	}
};


struct ElementBounds {
	double x_min, x_max, y_min, y_max, z_min, z_max;
	int original_index;
};