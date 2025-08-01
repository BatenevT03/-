//---------------------------------------------------------------------------
#include "myclass.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <omp.h> 
#include <array>
using namespace std;

bool parallel = false;

// [номер_параллелепипеда][точки_Гаусса_в_нём]
// Хранит значения базисных функций для всех элементов и точек Гаусса
// Формат: [номер элемента][номер точки Гаусса][номер базисной функции]
std::vector<std::vector<std::array<double, 8>>> phi_at_gauss;



vector<vector<double>> c = { {8},{4,8},{4,2,8},{2,4,4,8},{4,2,2,1,8},{2,4,1,2,4,8},{2,1,4,2,4,2,8},{1,2,2,4,2,4,4,8} };

vector<vector<double>> bx = { {4},{-4,4},{2,-2,4},{-2,2,-4,4},{2,-2,1,-1,4,4},{-2,2,-1,1,-4,4},{1,-1,2,-2,2,-2,4},{-1,1,-2,2,-2,2,-4,4} };

vector<vector<double>> by = { {4},{2,4},{-4,-2,4},{-2,-4,2,4},{2,1,-2,-1,4},{1,2,-1,-2,2,4},{-2,-1,2,1,-4,-2,4},{-1,-2,1,2,-2,-4,2,4} };

vector<vector<double>> bz = { {4},{2,4},{2,1,4},{1,2,2,4},{-4,-2,-2,-1,4},{-2,-4,-1,-2,2,4},{-2,-1,-4,-2,2,1,4},{-1,-2,-2,-4,1,2,2,4} };
vector<double> right_vector(M);
vector<vector<double>> c1 = { {4},{2,4},{2,1,4},{1,2,2,4} };
vector<vector<double>> rf_v = { {-4,4,-2,2,-2,2,-1,1},{-4,4,-2,2,-2,2,-1,1},{-2,2,-4,4,-1,1,-2,2},{-2,2,-4,4,-1,1,-2,2} };


// Константы квадратуры Гаусса
constexpr array<double, 3> tau = { 8. / 9., 5. / 9., 5. / 9. };
constexpr array<double, 3> t = { 0., 0.77459666924148337, -0.77459666924148337 };

double right(vector<double>& x)
{
	//if (x[0] == 0. && x[1] == 0. && x[2] == 0) return 1.0;
	//if (x[0] == 0. && x[1] == 0. && x[2] == 0) return 1000.;

	return 0.;
	return -6;
	return (-6 * x[0] - 6 * x[1] - 6 * x[2]);
	return sin(x[0] + x[1] + x[2]);
}

double dudn(vector<double>& x, int i)
{
	//return 1;
	//return cos(x[0] + x[1] + x[2]);

	switch (i)
	{
	case 1:
		return 2 * x[0];
		return 3 * x[0] * x[0];
	case 2:
		return 2 * x[1];
		return 3 * x[1] * x[1];
	case 3:
		return 2 * x[2];
		return 3 * x[2] * x[2];
	}

}



int u(int i)
{
	return (i ) % 2 + 1-1;
}
int v(int i)
{
	return ((i ) / 2) % 2 + 1-1;
}
int vu(int i)
{
	return (i ) / 4 + 1-1;
}




// Базисные функции и их производные
inline double dphidxi(double eta, double theta, int i) {
	const double mt = 1 - theta, t_ = theta;
	const double me = 1 - eta, e = eta;
	switch (i) {
	case 0: return -me * mt;
	case 1: return  me * mt;
	case 2: return -e * mt;
	case 3: return  e * mt;
	case 4: return -me * t_;
	case 5: return  me * t_;
	case 6: return -e * t_;
	case 7: return  e * t_;
	default: return 0;
	}
}

inline double dphideta(double xi, double theta, int i) {
	const double mt = 1 - theta, t_ = theta;
	const double mx = 1 - xi, x = xi;
	switch (i) {
	case 0: return -mx * mt;
	case 1: return -x * mt;
	case 2: return  mx * mt;
	case 3: return  x * mt;
	case 4: return -mx * t_;
	case 5: return -x * t_;
	case 6: return  mx * t_;
	case 7: return  x * t_;
	default: return 0;
	}
}

inline double dphidtheta(double xi, double eta, int i) {
	const double mx = 1 - xi, x = xi;
	const double me = 1 - eta, e = eta;
	switch (i) {
	case 0: return -mx * me;
	case 1: return -x * me;
	case 2: return -mx * e;
	case 3: return -x * e;
	case 4: return  mx * me;
	case 5: return  x * me;
	case 6: return  mx * e;
	case 7: return  x * e;
	default: return 0;
	}
}

struct Derivatives {
	double dxdxi, dxdeta, dxdtheta;
	double dydxi, dydeta, dydtheta;
	double dzdxi, dzdeta, dzdtheta;
};

Derivatives computeDerivatives(const vector<double>& x, const vector<double>& y, const vector<double>& z,
	double xi, double eta, double theta) {
	Derivatives d{};
	for (int i = 0; i < 8; ++i) {
		double dxi_i = dphidxi(eta, theta, i);
		double deta_i = dphideta(xi, theta, i);
		double dtheta_i = dphidtheta(xi, eta, i);

		d.dxdxi += x[i] * dxi_i;
		d.dxdeta += x[i] * deta_i;
		d.dxdtheta += x[i] * dtheta_i;

		d.dydxi += y[i] * dxi_i;
		d.dydeta += y[i] * deta_i;
		d.dydtheta += y[i] * dtheta_i;

		d.dzdxi += z[i] * dxi_i;
		d.dzdeta += z[i] * deta_i;
		d.dzdtheta += z[i] * dtheta_i;
	}
	return d;
}

inline double computeJ(const Derivatives& d) {
	return (d.dxdxi * (d.dydeta * d.dzdtheta - d.dydtheta * d.dzdeta) -
		d.dxdeta * (d.dydxi * d.dzdtheta - d.dydtheta * d.dzdxi) +
		d.dxdtheta * (d.dydxi * d.dzdeta - d.dydeta * d.dzdxi));
}

double f(const Derivatives& d, double xi, double eta, double theta, int i, int j, int c) {
	double phixi = dphidxi(eta, theta, i);
	double phieta = dphideta(xi, theta, i);
	double phitheta = dphidtheta(xi, eta, i);

	double phixi2 = dphidxi(eta, theta, j);
	double phieta2 = dphideta(xi, theta, j);
	double phitheta2 = dphidtheta(xi, eta, j);

	switch (c) {
	case 1: {
		double a = phixi * (d.dydeta * d.dzdtheta - d.dydtheta * d.dzdeta) +
			phieta * (d.dydtheta * d.dzdxi - d.dydxi * d.dzdtheta) +
			phitheta * (d.dydxi * d.dzdeta - d.dydeta * d.dzdxi);
		double b = phixi2 * (d.dydeta * d.dzdtheta - d.dydtheta * d.dzdeta) +
			phieta2 * (d.dydtheta * d.dzdxi - d.dydxi * d.dzdtheta) +
			phitheta2 * (d.dydxi * d.dzdeta - d.dydeta * d.dzdxi);
		return a * b;
	}
	case 2: {
		double a = d.dxdxi * (phieta * d.dzdtheta - phitheta * d.dzdeta) +
			d.dxdeta * (phitheta * d.dzdxi - phixi * d.dzdtheta) +
			d.dxdtheta * (phixi * d.dzdeta - phieta * d.dzdxi);
		double b = d.dxdxi * (phieta2 * d.dzdtheta - phitheta2 * d.dzdeta) +
			d.dxdeta * (phitheta2 * d.dzdxi - phixi2 * d.dzdtheta) +
			d.dxdtheta * (phixi2 * d.dzdeta - phieta2 * d.dzdxi);
		return a * b;
	}
	case 3: {
		double a = d.dxdxi * (d.dydeta * phitheta - d.dydtheta * phieta) +
			d.dxdeta * (d.dydtheta * phixi - d.dydxi * phitheta) +
			d.dxdtheta * (d.dydxi * phieta - d.dydeta * phixi);
		double b = d.dxdxi * (d.dydeta * phitheta2 - d.dydtheta * phieta2) +
			d.dxdeta * (d.dydtheta * phixi2 - d.dydxi * phitheta2) +
			d.dxdtheta * (d.dydxi * phieta2 - d.dydeta * phixi2);
		return a * b;
	}
	default: return 0;
	}
}










double Gauss2(const vector<double>& x, const vector<double>& y, const vector<double>& z,
	int i, int j, int c, double a, double b) {
	double result = 0;
	const double mid = (a + b) * 0.5;
	const double hHalf = (b - a) * 0.5;

	for (int gi = 0; gi < 3; ++gi) {
		double xi = mid + t[gi] * hHalf;
		for (int gj = 0; gj < 3; ++gj) {
			double eta = mid + t[gj] * hHalf;
			for (int gk = 0; gk < 3; ++gk) {
				double theta = mid + t[gk] * hHalf;
				Derivatives d = computeDerivatives(x, y, z, xi, eta, theta);
				double jac = computeJ(d);
				double jacSqInv = 1.0 / (jac * jac);
				result += tau[gi] * tau[gj] * tau[gk] * f(d, xi, eta, theta, i, j, c)/jac;
			}
		}
	}
	return (b - a) * (b - a) * (b - a) * result * 0.125;
}

inline void mul_c_vector(vector<double>&vec, vector<double>&result)
{
	for(int i=0;i<8;i++) {
		result[i]=0;

		for(int j=0;j<8;j++) {
			int i1,j1;
			if(i>j) {
				i1=i;
				j1=j;
			}
			else {
				i1=j;
				j1=i;
			}
			result[i]+=c[i1][j1]*vec[j];
		}
	}
}

void mul_c1_vector(vector<double>&vec, vector<double>&result)
{
	for(int i=0;i<4;i++) {
		result[i]=0;

		for(int j=0;j<4;j++) {
			int i1,j1;
			if(i>j) {
				i1=i;
				j1=j;
			}
			else {
				i1=j;
				j1=i;
			}
			result[i]+=c1[i1][j1]*vec[j];
		}
	}	
}

double scalar(vector<double>&vector1, vector<double>&vector2,int n)
{
	double f=0;
	for(int i=0;i<n;i++)
		f+=vector1[i]*vector2[i];
	return f;
}

double norma(vector<double>& vector,int n)
{
	double f=scalar(vector,vector,n);
	return sqrt(f);
}

void sum(vector<double>&vector1, vector<double>&vector2, vector<double>&result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=vector1[i]+vector2[i];
}
//-----------------------------
void mul(double a, vector<double>&vec, vector<double>&result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=a*vec[i];
}

void mul2(vector<double>&a, vector<double>&vec, vector<double>&result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=a[i]*vec[i];
}
//-----------------------------
void mov(vector<double>&vec, vector<double>&result,int n)
{
	for(int i=0;i<n;i++)
		result[i]=vec[i];
}
//-----------------------------
int local::count=0;
int cross::count=0;
//----------------------------------
inline void MATRIX::mul_matrix_vector(vector<double>&vec, vector<double>&result)
{
	if (parallel)
	{
		//iteration++;
		// Проверка размеров
		if (vec.size() != n || result.size() != n) {
			throw std::invalid_argument("Размеры вектора и матрицы не совместимы");
		}

		// Инициализация результата нулями
		std::fill(result.begin(), result.end(), 0.0);

		// Умножение диагональных элементов
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			result[i] = d[i] * vec[i];
		}

		// Умножение недиагональных элементов
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			for (int j = ig[i]; j < ig[i + 1]; j++) {
				int col_index = jg[j - 1] - 1;  // Предполагаем, что индексация начинается с 1
				double value = gg[j - 1];

#pragma omp atomic
				result[i] += value * vec[col_index];

#pragma omp atomic
				result[col_index] += value * vec[i];
			}
		}
	}
	else
	{
		int i;
		for(i=0;i<n;i++)
			result[i]=d[i]*vec[i];

		for(i=0;i<n;i++)
			for(int j=ig[i];j<ig[i+1];j++) {
				result[i]+=gg[j-1]*vec[jg[j-1]-1];
				result[jg[j-1]-1]+=gg[j-1]*vec[i];
			}
	}
	

}
//-------------------------------------------------------
void MATRIX::solution_x_l(vector<double>&f, vector<double>&x)
{
    int i;
	for(i=0;i<n;i++)
		x[i]=f[i];
	//SY=F
	for(i=0;i<n;i++){
		for(int j=ig[i]-1;j<ig[i+1]-1;j++)
			x[i]-=gg[j]*x[jg[j]-1];
		x[i]/=d[i];
	}
}

//------------------------------
void MATRIX::solution_x_u(vector<double>&f, vector<double>&x)
{
    int i;
	for(i=0;i<n;i++)
		x[i]=f[i];

	for(i=n-1;i>=0;i--) {
		x[i]/=d[i];
		for(int j=ig[i+1]-2;j>=ig[i]-1;j--)
			x[jg[j]-1]-=gg[j]*x[i];
	}

}
//----------------------------------
void MATRIX::solution_x(vector<double>&f, vector<double>&x)
{
    int i;
	for (i = 0;i < n;i++)
	{
		x[i] = f[i];
	}

	//SY=F 
	for(i=0;i<n;i++) {
		for(int j=ig[i]-1;j<ig[i+1]-1;j++)
			x[i]-=gg[j]*x[jg[j]-1];
		x[i]/=d[i];
	}

	//STX=Y
	for(i=n-1;i>=0;i--){
		x[i]/=d[i];
		for(int j=ig[i+1]-2;j>=ig[i]-1;j--)
			x[jg[j]-1]-=gg[j]*x[i];
	}
}


//----------------------------------
void MATRIX::SST(MATRIX *my)
{
    int i,j,k,l;

    my->n=n;   
	my->d.resize(n);
	my->ig.resize(n + 1);
	


    for(i=0;i<=n;i++)
		my->ig[i]=ig[i];

	my->jg.resize(ig[n]-1);
	my->gg.resize(ig[n] - 1);

    for(i=0;i<ig[n]-1;i++)
		my->jg[i]=jg[i];

	for(j=0;j<n;j++)
	{
	double sum=0;
		for(k=ig[j]-1;k<ig[j+1]-1;k++)
		sum+=my->gg[k]*my->gg[k];

	my->d[j]=sqrt(fabs(d[j]-sum));

	for(i=j+1;i<n;i++)
	{
		int number;
		int flag=1;
		for(l=ig[i]-1;l<ig[i+1]-1&&flag;l++)
		if(jg[l]==j+1) flag=0;

		number=l-1;
		if(flag) continue;

		sum=0;
		for(k=ig[i]-1;k<ig[i+1]-1&&jg[k]<=j;k++)
		{
			flag=1;
			for(l=ig[j]-1;l<ig[i+1]-1&&flag&&jg[l]<=jg[k];l++)
			if(jg[l]==jg[k]) flag=0;

			l--;

			if(!flag)
			sum+=my->gg[l]*my->gg[k];

		}

	my->gg[number]=(gg[number]-sum)/my->d[j];
	}


	}

}
MATRIX::~MATRIX()
{
	//delete [] d;
	//delete [] ig;
	//delete [] gg;
	//delete [] jg;
}
 
//Считывание нумерации
void GLOBAL_MATRIX::read_local()
{
	ifstream in("local.txt");
	in >> local::count;
	cout << "Finite elements: " << local::count << endl;
	matr = new local[local::count];

	for(int i=0;i<local::count;i++) {
		for (int j = 0;j < 8;j++)
		{
			matr[i].mas.resize(8);
			in >> matr[i].mas[j];
		}
		in >> matr[i].lambda;
	}

	
}
//----------------------------------------------------
//Считывание координат узлов
void GLOBAL_MATRIX::read_cross()
{
	//set == point
	ifstream in("cross.txt");
	in >> cross::count;
	n=cross::count;	
	set = new cross[n];
	cout << "Cross count: " << n << endl;
	for (int i = 0;i < cross::count;i++)
		for (int j = 0;j < 3;j++)
		{
			set[i].mas.resize(3);
			in >> set[i].mas[j];
		}
	/*for(int i=0;i<cross::count;i++)
	{
	printf("\n");
	for(int j=0;j<3;j++)
	printf("%2.1f ",set[i].mas[j]);
	} */
	in.close();
}




//void GLOBAL_MATRIX::formier_matrix()//параллелепипеды
//{
//	double h[3];
//	for (int i = 0; i < n; i++) {
//		d[i] = 0.0;
//		f[i] = 0.0;
//	}
//	for (int i = 0; i < ig[n] - 1; i++) {
//		gg[i] = 0.0;
//	}
//	for (int k = 0;k < local::count;k++)
//	{
//		//Вычисление шага
//		//h----------------------------------------------------------
//		for (int i = 0;i < 3;i++) { //i-по x,y,z    
//			int flag = 1;
//			int j;
//			for (j = 1; j < M && flag;j++) {//1 узел фиксируем,пробегаем по остальным    
//				flag = 0;
//				for (int l = 0;l < 3 && !flag;l++)//проверяем,лежат ли точки на нужном ребре
//					if (i != l && set[matr[k].mas[0] - 1].mas[l] != set[matr[k].mas[j] - 1].mas[l])
//						flag = 1;
//			}
//			if (!flag)
//				h[i] = fabs(set[matr[k].mas[0] - 1].mas[i] - set[matr[k].mas[j - 1] - 1].mas[i]);
//		}
//		//------------------------------------------------------------
//		//формирование элементов матрицы
//		//заполнение мссива gg
//		double b_k = matr[k].lambda * h[0] * h[1] * h[2] / 36.;
//		double c_k = h[0] * h[1] * h[2] / 216.;
//		//!с_k заменить
//
//		//!
//		//double c_k=1,b_k=0;
//		vector<double> fr(M);
//		//вектор правой части для локал. матрицы
//		for (int i = 0;i < M;i++)
//			right_vector[i] = right(set[matr[k].mas[i] - 1].mas);
//		mul_c_vector(right_vector, fr);
//
//		//cout << "Right vector in element " << k << ": ";
//		//for (int i = 0; i < M; i++) cout << right_vector[i] << " ";
//		//cout << endl;
//
//
//		for (int i = 0;i < M;i++) {
//			for (int j = 0;j < i;j++) {
//				int i1 = matr[k].mas[i] - 1;
//				int j1 = matr[k].mas[j] - 1;
//				double s = b_k / (h[0] * h[0]) * bx[i][j] +
//					b_k / (h[1] * h[1]) * by[i][j] + b_k / (h[2] * h[2]) * bz[i][j];
//
//				//добавка в gg
//				if (i1 < j1)
//					add(j1, i1, s);
//				else
//					add(i1, j1, s);
//			}
//			//добавка в диагональ
//			d[matr[k].mas[i] - 1] += b_k / (h[0] * h[0]) * bx[i][i] +
//				b_k / (h[1] * h[1]) * by[i][i] + b_k / (h[2] * h[2]) * bz[i][i];
//			//добавка в правую часть
//			f[matr[k].mas[i] - 1] += c_k * fr[i];
//			//!
//
//			//!
//		}
//	}
//	ofstream q("d.txt");
//	for (int i = 0;i < d.size();i++)
//	{
//		q << d[i] << endl;
//	}
//}


void GLOBAL_MATRIX::formier_matrix() {//шестигранники
	

	// Инициализация матрицы (ваш существующий код)
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		f[i] = 0.0;
	}
	for (int i = 0; i < ig[n] - 1; i++) {
		gg[i] = 0.0;
	}

	// Цикл по элементам
	for (int k = 0; k < local::count; k++) {
		//cout << "Is: " << k << endl;
		vector<double> x(8), y(8), z(8);
		for (int i = 0; i < 8; i++) {
			x[i] = set[matr[k].mas[i] - 1].mas[0];
			y[i] = set[matr[k].mas[i] - 1].mas[1];
			z[i] = set[matr[k].mas[i] - 1].mas[2];
		}

		// Обнуляем локальную матрицу
		vector<vector<double>> Ke(8, vector<double>(8, 0.0));

		// Цикл по предвычисленным точкам Гаусса
		

			// Вычисляем глобальные производные
			for (int i = 0; i < 8; ++i) {
				for (int j = 0; j < 8; ++j) {
					Ke[i][j] = matr[k].lambda*(Gauss2(x, y, z, i, j, 1, 0, 1) +
						Gauss2(x, y, z, i, j, 2, 0, 1) +
						Gauss2(x, y, z, i, j, 3, 0, 1));
					//cout << Ke[i][j] << " ";
				}
				//cout << endl;
			}
		

		// Добавление в глобальную матрицу (ваш существующий код)
		for (int i = 0; i < 8; i++) {
			int global_i = matr[k].mas[i] - 1;
			d[global_i] += Ke[i][i];

			for (int j = 0; j < i; j++) {
				int global_j = matr[k].mas[j] - 1;
				double value = Ke[i][j];

				if (fabs(value) < 1e-16) value = 0;

				if (global_i < global_j) {
					add(global_j, global_i, value);
				}
				else {
					add(global_i, global_j, value);
				}
			}
		}
	}
	/*ofstream q("d.txt");
	for (int i = 0;i < d.size();i++)
	{
		q << d[i] << endl;
	}*/
}


//------------------------------------------------------
void GLOBAL_MATRIX::add(int i,int j,double x)
{
	int k;
	for(k=ig[i]-1;k<ig[i+1]-1;k++)
		if(jg[k]==j+1)
			break;
	gg[k]+=x;
}


void GLOBAL_MATRIX::formier_profil() {
	std::vector<std::pair<int, int>> list; // Временный массив для хранения списка
	std::vector<int> listbeg(n + 1, 0); // Временный массив для хранения начал списков (индексация с 1)
	int listsize = 0; // Текущая длина массива list

	// Цикл формирования списка одним проходом по элементам
	for (int ielem = 0; ielem < local::count; ++ielem) {
		int numUnknowns = 8;
		for (int i = 0; i < numUnknowns; ++i) {
			int k = matr[ielem].mas[i];
			for (int j = i + 1; j < numUnknowns; ++j) {
				int ind1 = k;
				int ind2 = matr[ielem].mas[j];
				if (ind2 < ind1) {
					std::swap(ind1, ind2);
				}

				// Заносим связь большего номера с меньшим, т.е. ind2 с ind1
				int iaddr = listbeg[ind2];
				if (iaddr == 0) {
					// Списка ещё не было, создаём его
					listsize++;
					listbeg[ind2] = listsize;
					list.push_back({ ind1, 0 });
				}
				else {
					// Список уже был. Ищем в нём ind1
					while (list[iaddr - 1].first < ind1 && list[iaddr - 1].second > 0) {
						iaddr = list[iaddr - 1].second;
					}
					if (list[iaddr - 1].first > ind1) {
						listsize++;
						list.push_back({ list[iaddr - 1].first, list[iaddr - 1].second });
						list[iaddr - 1] = { ind1, listsize };
					}
					else if (list[iaddr - 1].first < ind1) {
						// Не нашли, а список закончился. Добавляем в конец списка
						listsize++;
						list[iaddr - 1].second = listsize;
						list.push_back({ ind1, 0 });
					}
				}
			}
		}
		
	}

	// Создание портрета по списку
	ig.resize(n + 1); // Размер N+1, чтобы индексация начиналась с 0
	ig[0] = 1; // ig[0] инициализируется единицей
	for (int i = 1; i <= n; ++i) {
		ig[i] = ig[i - 1];
		int iaddr = listbeg[i];
		while (iaddr != 0) {
			jg.push_back(list[iaddr - 1].first);
			ig[i]++;
			iaddr = list[iaddr - 1].second;
		}
	}
	//for (int i = 0; i < jg.size(); ++i) {
	//	printf("%d ", jg[i]);
	//}
	gg.resize(jg.size());
}





GLOBAL_MATRIX::GLOBAL_MATRIX()
{
	
	POR=1e30;
	//generate();
	read_cross();
	read_local();

	d.resize(n); 
	f.resize(n); 
	x.resize(n); 
	ig.resize(n + 1);

	formier_profil();
	formier_matrix();
	int flag = 0;
	cout << "cross cout: " << cross::count << endl;
	ifstream in("el.txt");
	vector<double> a(3, 0), b(3, 0), c(3, 0);
	in >> a[0] >> a[1] >> a[2];
	in >> b[0] >> b[1] >> b[2];
	in >> c[0] >> c[1] >> c[2];
	for (int i = 0;i < cross::count;i++)
	{
		if (set[i].mas[0] == a[0] && set[i].mas[1] == a[1] && set[i].mas[2] == a[2]) f[i] = 1.;//общий
		else if (set[i].mas[0] == b[0] && set[i].mas[1] == b[1] && set[i].mas[2] == b[2]) f[i] = 1.;//правый
		else if ((set[i].mas[0] == c[0] && set[i].mas[1] == c[1] && set[i].mas[2] == c[2])) f[i] = 1.;//нижний
		else f[i] = 0.;
	}
	in.close();
	ofstream out("f.txt");
	for (int i = 0;i < cross::count;i++)
	{
		out << f[i] << endl;
	}
	out.close();
	



	//for (int i = cross::count/2;i < cross::count;i++)
	//{
	//	
	//	if (set[i].mas[0] == 0. && set[i].mas[1] == 0. && set[i].mas[2] == 0.)
	//	{
	//		f[i] = 1000.;
	//		break;
	//	}
	//}
	//printf("\n=============\n");
	//for(int i=0;i<ig[n]-1;i++)
	//	printf("%1.2f ",gg[i]);
	//printf("\n=============\n");
	//for(int i=0;i<n;i++)
	//	printf("%1.2f ",f[i]);
	//printf("\n=============\n");
	//for(int i=0;i<n;i++)
	//	printf("%1.2f ",d[i]);
	//printf("\n=============\n");
} 

GLOBAL_MATRIX::~GLOBAL_MATRIX()
{
	delete [] matr;
	delete [] set;
	f.clear();
	x.clear();
}
//---------------------------------------------
void GLOBAL_MATRIX::MSG()
{

	double e=1e-25;//-70
	int max=1000;
	MATRIX s;
	SST(&s);
	
	vector<double> r  (n,0);
	vector<double> z  (n,0);
	vector<double> h  (n,0);
	vector<double> h2 (n,0);
	vector<double> h3 (n,0);

	for(int i=0;i<n;i++) 
		x[i]=0;
	mul_matrix_vector(x,h);//h=Ax
	mul(-1,h,h,n);
	sum(f,h,r,n);//r=f-Ax

	s.solution_x(r,z);//z=M^-1*r
	//cout << norma(r, n) << " " << norma(f, n);
	double e2=norma(r,n)/norma(f,n);

	double a,b;
	int k;
	for(k=1;k<max&&e2>e;k++) {
		//alpha
		s.solution_x(r,h2);//h2=M^-1*r_k-1
		mul_matrix_vector(z,h);//h=A*z_k-1
		a=scalar(h2,r,n)/scalar(h,z,n);

		//r
		mov(r,h3,n);
		mul(-a,h,h,n);
		sum(r,h,r,n);

		//x
		mul(a,z,h,n);
		sum(x,h,x,n);

		//beta
		b=1/scalar(h2,h3,n);
		s.solution_x(r,h2);
		b*=scalar(h2,r,n);

		//z
		mul(b,z,z,n);
		sum(z,h2,z,n);

		//for (int i = 0;i < n;i++)
		//	printf("%20.18e\n", x[i]);
		e2=norma(r,n)/norma(f,n);
	}
	printf("k=%d e=%10.8e\n",k-1,e2);	

	r.clear();
	z.clear();
	h.clear();
	h2.clear();
	h3.clear();
	//printf("\n\n\n%d\n\n\n", iteration);
}
//--------------------------------------

//---------------------------------------
void GLOBAL_MATRIX::kraev_1()
{
	int i;
	double q;
	ifstream in("kraev_1.txt");
	int flag;
	while (in >> i >> q)
	{
		i--;
		d[i]=POR;
		f[i]=q*POR;
	}
	in.close();
	//for (int i = 0; i < ig.size(); i++)
	//{
	//	ig[i]--;
	//}
	//for (int i = 0; i < jg.size(); i++)
	//{
	//	jg[i]--;
	//}
}


//----------------------------------------
void GLOBAL_MATRIX::kraev_2()
{
	vector<int> cr(4,0);//граничные узлы
	vector<double> teta(4,0);//значения потока в узлах
	ifstream in("kraev_2.txt");
	while(in>>cr[0]>>cr[1]>>cr[2]>>cr[3]>>teta[0]>>teta[1]>>teta[2]>>teta[3]) 
	{
		double hxyz = 1.0;
		double eps=1e-12;
		double hx = 0.0, hy = 0.0;
		for(int i = 0; i < 3; i++) {
			double tmp = fabs(set[cr[0]-1].mas[i]-set[cr[1]-1].mas[i]);
			if(tmp > eps)
				hx += tmp;
			tmp = fabs(set[cr[0]-1].mas[i]-set[cr[2]-1].mas[i]);
			if(tmp > eps)
				hy += tmp;
		}
		hxyz = hx * hy;

		vector<double> vec_tet(4,0);
		mul_c1_vector(teta,vec_tet);
		//cout << endl;
		for (int i = 0;i < 4;i++)
		{
			//cout <<"Do "<< f[cr[i] - 1] << " ";
			f[cr[i] - 1] += hxyz / 36.0 * vec_tet[i];
			//cout <<"Posle "<< f[cr[i] - 1] << " ";
		}
	}
	in.close();
}



double GLOBAL_MATRIX::get_psi(int num_fe, int num_basis, double x, double y, double z)
{
	// Получаем узлы элемента
	int n1 = matr[num_fe].mas[0];  // (x0,y0,z0) - "нижний-левый-ближний" узел
	int n2 = matr[num_fe].mas[1];  // (x1,y0,z0)
	int n3 = matr[num_fe].mas[2];  // (x0,y1,z0)
	int n5 = matr[num_fe].mas[4];  // (x0,y0,z1)

	// Вычисляем шаги сетки
	double hx = set[n2 - 1].mas[0] - set[n1 - 1].mas[0];
	double hy = set[n3 - 1].mas[1] - set[n1 - 1].mas[1];
	double hz = set[n5 - 1].mas[2] - set[n1 - 1].mas[2];

	// Нормированные координаты относительно первого узла
	double xi = (x - set[n1 - 1].mas[0]) / hx;
	double eta = (y - set[n1 - 1].mas[1]) / hy;
	double zeta = (z - set[n1 - 1].mas[2]) / hz;

	// Одномерные линейные базисные функции
	double X1 = 1.0 - xi;
	double X2 = xi;
	double Y1 = 1.0 - eta;
	double Y2 = eta;
	double Z1 = 1.0 - zeta;
	double Z2 = zeta;

	// Трилинейные базисные функции для каждого узла
	switch (num_basis) {
	case 0: return X1 * Y1 * Z1;  // Узел n1 (x0,y0,z0)
	case 1: return X2 * Y1 * Z1;  // Узел n2 (x1,y0,z0)
	case 2: return X1 * Y2 * Z1;  // Узел n3 (x0,y1,z0)
	case 3: return X2 * Y2 * Z1;  // Узел n4 (x1,y1,z0)
	case 4: return X1 * Y1 * Z2;  // Узел n5 (x0,y0,z1)
	case 5: return X2 * Y1 * Z2;  // Узел n6 (x1,y0,z1)
	case 6: return X1 * Y2 * Z2;  // Узел n7 (x0,y1,z1)
	case 7: return X2 * Y2 * Z2;  // Узел n8 (x1,y1,z1)
	default: return 0.0;
	}
}
//----------------------------------------
double GLOBAL_MATRIX::U(double nx, double y, double z)
{
	bool flag = false;
	int i;
	//находим элемент
	for(i = 0; i < local::count && !flag; i++) {
		int v1  = matr[i].mas[0];
		int v2  = matr[i].mas[1];
		int v3  = matr[i].mas[2];
		int v5  = matr[i].mas[4];
		//проверка, чтобы в i-ом КЭ наша точка лежала внутри него
		if(nx >= set[v1-1].mas[0] && nx <= set[v2-1].mas[0] &&
			y >= set[v1-1].mas[1] && y <= set[v3-1].mas[1] &&
			z >= set[v1-1].mas[2] && z <= set[v5-1].mas[2])
			flag = true;
	}
	int lc = i;
	double sum = 0;
	//f(x0)=sum(q_i*psi_i(x0))
	for(i = 0; i < 8; i++)
		sum += x[matr[lc-1].mas[i]-1] * get_psi(lc-1,i,nx,y,z);
	return sum;
}





void GLOBAL_MATRIX::generate()
{
	int hx = 0, hy = 0, hz = 0;
	double x0 = 0, x1 = 1, h = (double)x1/2;
	ofstream in("cross.txt");
	//cin >> x0 >> x1 >> h;
	double l = pow(((x1 - x0) / h), 3);//число КЭ
	double n = pow(((x1 - x0) / h) + 1, 3);//число узлов
	double nx = (x1 - x0) / h + 1;
	double lx = (x1 - x0) / h;
	int i, j, k;
	in << n << endl << endl;
	for (i = 0;i < nx;i++)
	{
		for (j = 0;j < nx;j++)
		{
			for (k = 0;k < nx;k++)
			{
				in << k * h << " " << j * h << " " << i * h << endl;
			}
		}
	}
	in.close();

	double numx = 0, numy = 0, numz = 0;
	int temp = 1;
	in.open("local.txt");
	in << l << endl;
	for (i = 0;i < lx;i++, temp = 1 + i * nx * nx)
	{
		for (j = 0;j < lx;j++, temp++)
		{
			for (k = 0;k < lx;k++, temp++)
			{
				in << temp << " " << temp + 1 << " " << temp + nx << " " << temp + nx + 1 << " " << temp + nx * nx << " " << temp + nx * nx + 1 << " " << temp + nx * nx + nx << " " << temp + nx * nx + nx + 1 << endl;
				in << 1 << endl;
			}
		}
	}

	in.close();

}




