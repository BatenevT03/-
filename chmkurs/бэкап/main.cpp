#include "myclass.h"
#include "CGM.h"
#include <chrono>
using namespace std;

double ug(vector<double> x)
{
	//return 0.;
	return x[0] + x[1] + x[2];
	//return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
	//return sin(x[0] + x[1] + x[2]);
}


void GLOBAL_MATRIX::print_result()
{
	double* x_ = new double[n], o[3];
	int nx = pow(n, (double)1 / 3), count = 0;
	int i, j, k;
	double h = 1;
	double a, b, c, temp, error = 0.;
	ifstream in("cross.txt");
	in >> a;
	std::cout << a;
	FILE* f = fopen("out.txt", "w");
	for (int i = 0;i < n;i++)
	{
		in >> a >> b >> c;
		temp = fabs(x[i] - ug(set[i].mas));
		error += temp * temp;
		fprintf(f, "%.4e %.4e %.4e %.4e \n", a, b, c, x[i]);
	}
	printf("error: %e\n", sqrt(error));
	fclose(f);
}

void print(vector<double> vec)
{
	for (int i = 0;i < vec.size();i++)
	{
		printf("%e ", vec[i]);
	}
	printf("\n");
}

void print(vector<int> vec)
{
	for (int i = 0;i < vec.size();i++)
	{
		printf("%d ", vec[i]);
	}
	printf("\n");
}

vector<double> sorted(vector<double> vec)
{
	vector<double> t;
	int a = vec.size() - 1, b = 0;
	for (int i = 0;i < vec.size();i++)
	{
		t.push_back(vec[a]);
		a--;
	}
	return t;
}

void generate(double start, double end, double h0, double ratio)
{
	ofstream in("cross.txt");
	in << std::fixed << setprecision(16) << endl;
	int iter = 2, i, j, k;
	vector<double>h, grid;
	h.push_back(h0);
	grid.push_back(start);
	grid.push_back(start = start + h0);
	double t = h0;
	while (true)
	{
		t *= ratio;
		start += t;
		if (start > end) break;
		h.push_back(t);
		grid.push_back(start);
		iter++;
	}
	print(h);
	//print(grid);
	//printf("%d\n", iter);

	//printf("In file cross.txt:\n");
	vector<double> inverted;
	for (i = 0;i < grid.size();i++)
		inverted.push_back(-grid[i]);
	inverted = sorted(inverted);
	inverted.pop_back();
	vector<double>Grid;
	int Iter = 2 * iter - 1;
	in << Iter * Iter * Iter << endl;;
	//printf("%d\n", Iter*Iter*Iter);
	Grid.reserve(2 * inverted.size()); // ќптимизаци€: резервируем пам€ть заранее
	copy(inverted.begin(), inverted.end(), back_inserter(Grid));
	copy(grid.begin(), grid.end(), back_inserter(Grid));
	print(Grid);
	vector<vector<double>> G;
	int count = 0;
	for (i = 0;i < Iter;i++)
	{
		for (j = 0;j < Iter;j++)
		{
			for (k = 0;k < Iter;k++, count++)
			{
				in << Grid[k] << " " << Grid[j] << " " << Grid[i] << endl;
				//printf("%f %f %f\n", Grid[k], Grid[j], Grid[i]);
				G.push_back({ Grid[k], Grid[j], Grid[i] });
			}
		}
	}

	in.close();
	in.open("local.txt");
	//printf("In file local.txt:\n");

	int l = (Iter - 1) * (Iter - 1) * (Iter - 1), lx = Iter - 1;
	int temp = 1;
	in << l << endl;
	//printf("%d\n",l);//in << l << endl;
	//lx=iter
	for (i = 0;i < lx;i++, temp = 1 + i * Iter * Iter)
	{
		for (j = 0;j < lx;j++, temp++)
		{
			for (k = 0;k < lx;k++, temp++)
			{
				in << temp << " " << temp + 1 << " " << temp + Iter << " " << temp + Iter + 1 << " " << temp + Iter * Iter << " " << temp + Iter * Iter + 1 << " " << temp + Iter * Iter + Iter << " " << temp + Iter * Iter + Iter + 1 << endl;
				in << 1 << endl;
				//printf("%d %d %d %d %d %d %d %d\n", temp, temp + 1, temp + Iter, temp + Iter + 1, temp + Iter * Iter, temp + Iter * Iter + 1, temp + Iter * Iter + Iter, temp + Iter * Iter + Iter + 1);
				//printf("1\n");//in << 1 << endl;
			}
		}
	}
	in.close();















	in.open("kraev_1.txt");
	int a, b, c, d;
	//нижн€€ грань
	temp = 1;
	lx = Iter - 1;
	for (i = 0;i < lx;i++, temp++)
		for (j = 0;j < lx;j++, temp++)
		{
			a = temp;
			b = temp + 1;
			c = temp + Iter;
			d = temp + Iter + 1;
			//достаЄм точки через .mas, a,b,c,d номера 
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in << endl;
	//лева€ грань
	temp = 1;
	for (i = 0;i < lx;i++, temp += Iter)
		for (j = 0;j < lx;j++, temp += Iter)
		{
			a = temp;
			b = temp + Iter;
			c = temp + Iter * Iter;
			d = temp + Iter * Iter + Iter;
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in << endl;
	//передн€€ грань
	temp = 1;
	for (i = 0;i < lx;i++, temp = 1 + i * Iter * Iter)
		for (j = 0;j < lx;j++, temp++)
		{
			a = temp;
			b = temp + 1;
			c = temp + Iter * Iter;
			d = c + 1;
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in << endl;
	//верхн€€
	temp = 1 + Iter * Iter * lx;
	for (i = 0;i < lx;i++, temp++)
		for (j = 0;j < lx;j++, temp++)
		{
			a = temp;
			b = temp + 1;
			c = temp + Iter;
			d = c + 1;
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;

		}
	in << endl;
	//права€ грань
	temp = Iter;
	for (i = 0;i < lx;i++, temp += Iter)
		for (j = 0;j < lx;j++, temp += Iter)
		{
			a = temp;
			b = temp + Iter;
			c = temp + Iter * Iter;
			d = temp + Iter * Iter + Iter;
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in << endl;
	//задн€€ грань
	temp = 1 + lx * Iter;
	for (i = 0;i < lx;i++, temp = 1 + lx * Iter + Iter * Iter * i)
		for (j = 0;j < lx;j++, temp++)
		{
			a = temp;
			b = temp + 1;
			c = temp + Iter * Iter;
			d = c + 1;
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in.close();
}





int main()
{
	double begin = 0;
	double end = 1;
	double h0 = 0.3;
	double ratio = 1.0; //  оэффициент разр€дки (1.1 - уплотнение у начала)
	vector<double> abcd;
	//generate(begin, end, h0, ratio);
	cout << endl;
	//generate();
	GLOBAL_MATRIX my;
	//.generate();
	//print(my.ig);
	//my.kraev_2();
	my.kraev_1();
	
	//CGM r;
	//r.init(my.ig, my.jg, my.d, my.gg, my.f, my.n);
	//r.solve(my.x);
	//auto start_serial = chrono::high_resolution_clock::now();
	my.MSG();
	//auto end_serial = chrono::high_resolution_clock::now();
	//auto duration_serial = chrono::duration_cast<chrono::milliseconds>(end_serial - start_serial);
	// ¬ывод времени выполнени€
	//cout << "Time: " << duration_serial.count() << " milliseconds" << endl;
	my.print_result();
	double res;
	res = my.U(0.5, 0.5, 0.5);
	cout << endl << res << endl;
	return 0;
}