//---------------------------------------------------------------------------
#include "myclass.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <omp.h> 

using namespace std;

bool parallel = false;


vector<vector<double>> c = { {8},{4,8},{4,2,8},{2,4,4,8},{4,2,2,1,8},{2,4,1,2,4,8},{2,1,4,2,4,2,8},{1,2,2,4,2,4,4,8} };

vector<vector<double>> bx = { {4},{-4,4},{2,-2,4},{-2,2,-4,4},{2,-2,1,-1,4,4},{-2,2,-1,1,-4,4},{1,-1,2,-2,2,-2,4},{-1,1,-2,2,-2,2,-4,4} };

vector<vector<double>> by = { {4},{2,4},{-4,-2,4},{-2,-4,2,4},{2,1,-2,-1,4},{1,2,-1,-2,2,4},{-2,-1,2,1,-4,-2,4},{-1,-2,1,2,-2,-4,2,4} };

vector<vector<double>> bz = { {4},{2,4},{2,1,4},{1,2,2,4},{-4,-2,-2,-1,4},{-2,-4,-1,-2,2,4},{-2,-1,-4,-2,2,1,4},{-1,-2,-2,-4,1,2,2,4} };
vector<double> right_vector(M);
vector<vector<double>> c1 = { {4},{2,4},{2,1,4},{1,2,2,4} };
vector<vector<double>> rf_v = { {-4,4,-2,2,-2,2,-1,1},{-4,4,-2,2,-2,2,-1,1},{-2,2,-4,4,-1,1,-2,2},{-2,2,-4,4,-1,1,-2,2} };


double tau[3] = { 8. / 9.,5. / 9.,5. / 9. };
double t[3] = { 0.,0.77459666924148337,-0.77459666924148337 };

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

double dphidxi(double eta, double theta, int i)
{
	switch (i + 1)
	{
	case 1: return -(1 - eta) * (1 - theta);
	case 2: return (1 - eta) * (1 - theta);
	case 3: return -eta * (1 - theta);
	case 4: return eta * (1 - theta);
	case 5: return -(1 - eta) * theta;
	case 6: return (1 - eta) * theta;
	case 7: return -eta * theta;
	case 8: return eta * theta;
	default: return 0;
	}
}

double dphideta(double xi, double theta, int i)
{
	switch (i + 1)
	{
	case 1: return -(1 - xi) * (1 - theta);
	case 2: return -xi * (1 - theta);
	case 3: return (1 - xi) * (1 - theta);
	case 4: return xi * (1 - theta);
	case 5: return -(1 - xi) * theta;
	case 6: return -xi * theta;
	case 7: return (1 - xi) * theta;
	case 8: return xi * theta;
	default: return 0;
	}
}

double dphidtheta(double xi, double eta, int i)
{
	switch (i + 1)
	{
	case 1: return -(1 - xi) * (1 - eta);
	case 2: return -xi * (1 - eta);
	case 3: return -(1 - xi) * eta;
	case 4: return -xi * eta;
	case 5: return (1 - xi) * (1 - eta);
	case 6: return xi * (1 - eta);
	case 7: return (1 - xi) * eta;
	case 8: return xi * eta;
	default: return 0;
	}
}




double dxi(vector<double> a, double xi, double eta, double theta)//dy/dxi dx/dxi dz/dxi 
{
	double result = 0.;
	for (int i = 0;i < 8;i++)
	{
		result += a[i] * dphidxi(eta, theta, i);
	}
	return result;
}

double deta(vector<double> a, double xi, double eta, double theta)
{
	double result = 0.;
	for (int i = 0;i < 8;i++)
	{
		result += a[i] * dphideta(xi, theta, i);
	}
	return result;
}

double dtheta(vector<double> a, double xi, double eta, double theta)//dy/dxi dx/dxi dz/dxi 
{
	double result = 0.;
	for (int i = 0;i < 8;i++)
	{
		result += a[i] * dphidtheta(xi, eta, i);
	}
	return result;
}





double J(vector<double> x, vector<double> y, vector<double> z, double xi, double eta, double theta)
{

	double positive = dxi(x, xi, eta, theta) * deta(y, xi, eta, theta) * dtheta(z, xi, eta, theta) +
		deta(x, xi, eta, theta) * dtheta(y, xi, eta, theta) * dxi(z, xi, eta, theta) +
		dtheta(x, xi, eta, theta) * dxi(y, xi, eta, theta) * deta(z, xi, eta, theta);
	double negative = dtheta(x, xi, eta, theta) * deta(y, xi, eta, theta) * dxi(z, xi, eta, theta) +
		deta(x, xi, eta, theta) * dxi(y, xi, eta, theta) * dtheta(z, xi, eta, theta) +
		dxi(x, xi, eta, theta) * dtheta(y, xi, eta, theta) * deta(z, xi, eta, theta);
	return positive - negative;
}


double dphidx(vector<double> x, vector<double> y, vector<double> z, double xi, double eta, double theta, int i)
{


	return (dphidxi(eta, theta, i) * deta(y, xi, eta, theta) * dtheta(z, xi, eta, theta) +
		dphideta(xi, theta, i) * dtheta(y, xi, eta, theta) * dxi(z, xi, eta, theta) +
		dphidtheta(xi, eta, i) * dxi(y, xi, eta, theta) * deta(z, xi, eta, theta) -
		dphidtheta(xi, eta, i) * deta(y, xi, eta, theta) * dxi(z, xi, eta, theta) -
		dphideta(xi, theta, i) * dxi(y, xi, eta, theta) * dtheta(z, xi, eta, theta) -
		dphidxi(eta, theta, i) * dtheta(y, xi, eta, theta) * deta(z, xi, eta, theta)) / J(x, y, z, xi, eta, theta);
}

double dphidy(vector<double> x, vector<double> y, vector<double> z, double xi, double eta, double theta, int i)
{
	return (dxi(x, xi, eta, theta) * dphideta(xi, theta, i) * dtheta(z, xi, eta, theta) +
		deta(x, xi, eta, theta) * dphidtheta(xi, eta, i) * dxi(z, xi, eta, theta) +
		dtheta(x, xi, eta, theta) * dphidxi(eta, theta, i) * deta(z, xi, eta, theta) -
		dtheta(x, xi, eta, theta) * dphideta(xi, theta, i) * dxi(z, xi, eta, theta) -
		dxi(x, xi, eta, theta) * dphidtheta(xi, eta, i) * deta(z, xi, eta, theta) -
		deta(x, xi, eta, theta) * dphidxi(eta, theta, i) * dtheta(z, xi, eta, theta)) / J(x, y, z, xi, eta, theta);
}

double dphidz(vector<double> x, vector<double> y, vector<double> z, double xi, double eta, double theta, int i)
{
	return (dxi(x, xi, eta, theta) * deta(y, xi, eta, theta) * dphidtheta(xi, eta, i) +
		deta(x, xi, eta, theta) * dtheta(y, xi, eta, theta) * dphidxi(eta, theta, i) +
		dtheta(x, xi, eta, theta) * dxi(y, xi, eta, theta) * dphideta(xi, theta, i) -
		dtheta(x, xi, eta, theta) * deta(y, xi, eta, theta) * dphidxi(eta, theta, i) -
		deta(x, xi, eta, theta) * dxi(y, xi, eta, theta) * dphidtheta(xi, eta, i) -
		dxi(x, xi, eta, theta) * dtheta(y, xi, eta, theta) * dphideta(xi, theta, i)) / J(x, y, z, xi, eta, theta);
}

double f(vector<double> x, vector<double> y, vector<double> z, double xi, double eta, double theta, int i, int j, int c)
{
	switch (c)
	{
	case 1:
		return dphidx(x, y, z, xi, eta, theta, i) * dphidx(x, y, z, xi, eta, theta, j);
	case 2:
		return dphidy(x, y, z, xi, eta, theta, i) * dphidy(x, y, z, xi, eta, theta, j);
	case 3:
		return dphidz(x, y, z, xi, eta, theta, i) * dphidz(x, y, z, xi, eta, theta, j);
	default: return 0;
	}
	return 0;
}





double Gauss2(vector<double> x, vector<double> y, vector<double> z, int i1, int j1, int c, double a, double b)
{
	double xi, yi, zi, result = 0.;
	double h = b - a;
	for (int i = 0;i < 3;i++)
	{
		xi = (a + b) / 2 + t[i] * h / 2;
		for (int j = 0;j < 3;j++)
		{
			yi = (a + b) / 2 + t[j] * h / 2;
			for (int k = 0;k < 3;k++)
			{
				zi = (a + b) / 2 + t[k] * h / 2;
				result += tau[i] * tau[j] * tau[k] * f(x, y, z, xi, yi, zi, i1, j1, c);
			}
		}
	}
	return (b - a) * (b - a) * (b - a) * result / 8;
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






void GLOBAL_MATRIX::formier_matrix()
{
	double h[3];
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		f[i] = 0.0;
	}
	for (int i = 0; i < ig[n] - 1; i++) {
		gg[i] = 0.0;
	}
	for (int k = 0;k < local::count;k++)
	{
		vector<double>x(8, 0), y(8, 0), z(8, 0);
		for (int i = 0;i < M;i++)
		{
			x[i] = set[matr[k].mas[i] - 1].mas[0];
			y[i] = set[matr[k].mas[i] - 1].mas[1];
			z[i] = set[matr[k].mas[i] - 1].mas[2];
		}

		for (int i = 0;i < M;i++) {
			d[matr[k].mas[i] - 1] = (Gauss2(x, y, z, i, i, 1, 0, 1) + Gauss2(x, y, z, i, i, 2, 0, 1) + Gauss2(x, y, z, i, i, 3, 0, 1));
			for (int j = 0;j < i;j++) {
				int i1 = matr[k].mas[i] - 1;
				int j1 = matr[k].mas[j] - 1;
				double value = (Gauss2(x, y, z, i, j, 1, 0, 1) + Gauss2(x, y, z, i, j, 2, 0, 1) + Gauss2(x, y, z, i, j, 3, 0, 1));

				// Обнуление малых значений
				if (fabs(value) < 1e-16) 
					value = 0;
				//добавка в gg
				if (i1 < j1)
					add(j1, i1, value);
				else
					add(i1, j1, value);


			}

		}
	}
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
		if (set[i].mas[0] == a[0] && set[i].mas[1] == a[1] && set[i].mas[2] == a[2]) f[i] = -1.;//общий
		else if (set[i].mas[0] == b[0] && set[i].mas[1] == b[1] && set[i].mas[2] == b[2]) f[i] = 1.;//правый
		else if ((set[i].mas[0] == c[0] && set[i].mas[1] == c[1] && set[i].mas[2] == c[2])) f[i] = 0.;//нижний
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
	cout << norma(r, n) << " " << norma(f, n);
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




