#include "myclass.h"
#include "CGM.h"
#include <chrono>
#include <set>
#include <map>
using namespace std;

double ug(vector<double> x)
{
	//return 0.;
	//return x[0] + x[1] + x[2];
	//return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	return pow(x[0], 3) + pow(x[1], 3) + pow(x[2], 3);
	//return sin(x[0] + x[1] + x[2]);
}

void readgrid(string f1, string f2, string f3, vector<double>& gridx, vector<double>& gridy, vector<double>& gridz)
{
	ifstream in(f1);

	double t;
	for (int i = 0;!in.eof();i++)
	{
		in >> t;
		gridx.push_back(t);
	}
	in.close();
	in.open(f2);


	for (int i = 0;!in.eof();i++)
	{
		in >> t;
		gridy.push_back(t);
	}
	in.close();
	in.open(f3);


	for (int i = 0;!in.eof();i++)
	{
		in >> t;
		gridz.push_back(t);
	}
	in.close();
}


void refineGrid(vector<double>& grid, int split) {
	if (split == 0) {
		return; // Ничего не делаем
	}
	else if (split > 0) {
		// Уточнение сетки: добавляем split точек в каждом интервале
		vector<double> refined;
		for (size_t i = 0; i < grid.size() - 1; ++i) {
			refined.push_back(grid[i]);
			// Добавляем split точек между grid[i] и grid[i+1]
			for (int j = 1; j <= split; ++j) {
				double point = grid[i] + j * (grid[i + 1] - grid[i]) / (split + 1);
				refined.push_back(point);
			}
		}
		refined.push_back(grid.back());
		grid = refined;
	}
	else if (split == -1) {
		// Разрежение сетки: удаляем каждый второй узел
		if (grid.size() <= 2) return; // Нельзя разрежать дальше
		vector<double> coarsened;
		for (size_t i = 0; i < grid.size(); i += 2) {
			coarsened.push_back(grid[i]);
		}
		// Если было нечётное число узлов, добавляем последний
		//if (grid.size() % 2 != 0) {
			coarsened.push_back(grid.back());
		//}
			if (coarsened.back() == coarsened[coarsened.size() - 2]) 
				coarsened.pop_back();
		grid = coarsened;
	}
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

	//for (int i=0;i<11;i++)
	//	for (int j=0;j<11;j++)
	//		for (int k=0;k<11;k++)
	//			fprintf(f, "%.4e %.4e %.4e %.16e \n", (double)k, (double)j, (double)i, U6((double)k, (double)j, (double)i));

	for (int i = 0;i < n;i++)
	{
		in >> a >> b >> c;
		temp = fabs(x[i] - ug(set[i].mas));
		error += temp * temp;
		fprintf(f, "%.4e %.4e %.4e %.16e \n", a, b, c, x[i]);
	}
	fclose(f);
}



void GLOBAL_MATRIX::Receiver(ofstream& out, double Len, double start_x, double start_z)
{
	vector<double> gridx, gridy, gridz;
	int ind;
	ifstream in("cross.txt");
	int n, i, j;
	bool flag = true;
	in >> n;
	double t1, t2, t3;
	for (i = 0; i < n && flag; i++)
	{
		in >> t1 >> t2 >> t3;
		if (t2 == 0. && t3 == 1000.) flag = false;
	}
	for (j = i - 1; i < n; j++)
	{
		if (set[j].mas[1] != 0.) break;
	}
	for (int q = i - 1; q < j; q++)
	{
		gridx.push_back(set[q].mas[0]);
		gridy.push_back(set[q].mas[1]);
		gridz.push_back(set[q].mas[2]);
	}

	out << start_x << " " << 0.0 << " " << start_z << " " << U6(start_x, 0, start_z) << endl;
	cout << "From receiver:\n";

	const double target_step = 25.0;
	double total_distance = 0.0;
	double currx = start_x;
	double currz = start_z;
	double dx, dz, r, len = 0.;
	int I;
	for (I = 0; I < gridx.size(); I++)
	{
		if (gridx[I] > start_x)
		{
			i = I;
			break;
		}
	}
	cout << gridx[i] << " " << gridz[i] << endl;

	while (currx <= start_x + Len && i < gridx.size())
	{
		dx = gridx[i] - currx;
		dz = gridz[i] - currz;
		cout << dx << " " << dz << endl;
		if (abs(dz) < 1e-14) dz = 0.;
		r = sqrt(dx * dx + dz * dz);

		if (r >= target_step - len)
		{
			double ratio = (target_step - len) / r;
			currx += dx * ratio;
			currz += dz * ratio;

			out << currx << " " << 0.0 << " " << currz << " " << U6(start_x, 0, start_z) << endl;
			len = 0.;

			// Не уменьшаем i, так как мы могли не достичь текущей точки сетки
		}
		else
		{
			len += r;
			currx = gridx[i];
			currz = gridz[i];
			i++;
		}
	}

	cout << "End receive.\n";
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
	Grid.reserve(2 * inverted.size()); // Оптимизация: резервируем память заранее
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
	//нижняя грань
	temp = 1;
	lx = Iter - 1;
	for (i = 0;i < lx;i++, temp++)
		for (j = 0;j < lx;j++, temp++)
		{
			a = temp;
			b = temp + 1;
			c = temp + Iter;
			d = temp + Iter + 1;
			//достаём точки через .mas, a,b,c,d номера 
			in << a << " " << ug(G[a - 1]) << endl;
			in << b << " " << ug(G[b - 1]) << endl;
			in << c << " " << ug(G[c - 1]) << endl;
			in << d << " " << ug(G[d - 1]) << endl;
		}
	in << endl;
	//левая грань
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
	//передняя грань
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
	//верхняя
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
	//правая грань
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
	//задняя грань
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


void gen_file(string cross,string Set,vector<double>& gridx, vector<double>& gridy, vector<double>& gridz,bool flag)
{
	ofstream in(cross);
	in << std::fixed << setprecision(16) << endl;
	int nx = gridx.size(), ny = gridy.size(), nz = gridz.size(), N = nx * ny * nz;
	int i, j, k;
	in << N << endl;
	int count = 0;
	vector<vector<double>> G;
	for (i = 0;i < nz;i++)
		for (j = 0;j < ny;j++)
			for (k = 0;k < nx;k++)
			{
				G.push_back({ gridx[k], gridy[j], gridz[i] });
				in << gridx[k] << " " << gridy[j] << " " << gridz[i] << "\n";
			}

	in.close();
	in.open(Set);
	int temp = 1;
	in << (nx - 1) * (ny - 1) * (nz - 1) << endl;
	//printf("%d\n",l);//in << l << endl;
	//lx=iter
	for (i = 0;i < nz - 1;i++, temp = 1 + i * nx * ny)
		for (j = 0;j < ny - 1;j++, temp++)
			for (k = 0;k < nx - 1;k++, temp++)
			{
				in << temp << " " << temp + 1 << " " << temp + nx << " " << temp + nx + 1 << " " << temp + nx * ny << " " << temp + nx * ny + 1 << " " << temp + nx * ny + nx << " " << temp + nx * ny + nx + 1 << endl;
				//if (G[temp - 1][2] < 0.) in << 0.001 << endl;
				//if (0 <= G[temp - 1][2] && G[temp - 1][2] < 500) in << 0.02625 << endl;
				//if (500 <= G[temp - 1][2] && G[temp - 1][2] < 800) in << 0.07475 << endl;
				//if (800 <= G[temp - 1][2]) in << 0.1 << endl;
				in << 1 << endl;
			}

	in.close();
	if (flag)
	{
		in.open("kraev_1.txt");
		int a, b, c, d;



		//нижн¤¤ грань 1
		temp = 1;
		for (i = 0;i < nx * ny;i++, temp++)
		{
			in << temp << " " << ug(G[temp - 1]) << endl;
		}
		in << endl;
		//лева¤ грань 2
		temp = 1;
		for (i = 0;i < ny * nz;i++)
		{
			in << temp << " " << ug(G[temp - 1]) << endl;
			temp += nx;
		}
		in << endl;
		//передн¤¤ грань 3
		int g = 1;
		temp = 1;
		for (i = 0;i < nx * nz;i++)
		{
			in << temp << " " << ug(G[temp - 1]) << endl;
			if (g == nx)
			{
				temp += nx * ny - nx + 1;
				g = 1;
			}
			else
			{
				temp++;
				g++;
			}
		}
		in << endl;
		//верхн¤¤ 4
		temp = 1 + nx * ny * (nz - 1);
		for (i = 0;i < nx * ny;i++)
		{
		    in << temp << " " << ug(G[temp - 1]) << endl;
		    temp++;
		}
		in << endl;
		//права¤ грань
		temp = nx;
		for (i = 0;i < ny * nz;i++)
		{
			in << temp << " " << ug(G[temp - 1]) << endl;
			temp += nx;
		}
		in << endl;
		//задн¤¤ грань
		temp = nx * ny - nx + 1;
		g = 1;
		for (i = 0;i < nx * nz;i++)
		{
			in << temp << " " << ug(G[temp - 1]) << endl;
			if (g == nx) {
				temp += nx * ny - nx + 1;
				g = 1;
			}
			else
			{
				temp++;
				g++;
			}
		}
		in.close();
	}
}

void GLOBAL_MATRIX::print_receiver()
{
	ifstream in("receiver.txt");
	ofstream out("value_from_receiver.txt");
	out << fixed << setprecision(14);
	vector<vector<double>> arr;
	int nf1 = find6el(-43.47779299732488, 0, 1010.80644372056952);
	cout << "Reshenie v electrode 1:\n";
	for (int i = 0;i < 8 && nf1!=-1;i++)
	{
		double x = set[matr[nf1].mas[i]-1].mas[0];
		double y = set[matr[nf1].mas[i]-1].mas[1];
		double z = set[matr[nf1].mas[i]-1].mas[2];
		cout << x << " " << y << " " << z << " " << U6(x, y, z)<<endl;
	}
	cout << endl;
	nf1 = find6el(0, 0, 1012.0);
	cout << "Reshenie v electrode 2:\n";
	for (int i = 0;i < 8 && nf1!=-1;i++)
	{
		double x = set[matr[nf1].mas[i] - 1].mas[0];
		double y = set[matr[nf1].mas[i] - 1].mas[1];
		double z = set[matr[nf1].mas[i] - 1].mas[2];
		cout << x << " " << y << " " << z << " " << U6(x, y, z) << endl;
	}
	cout << endl;
	while (!in.eof())
	{
		double t1, t2, t3;
		in >> t1 >> t2 >> t3;
		//arr.push_back({ -t1,t2,t3 });
		out << t1 << " " << t3<<" " << U6(t1, t2, t3)<<endl;
	}
	//reverse(arr.begin(), arr.end());
	//for (int i = 0;i < arr.size();i++)
	//{
	//	out <<  arr[i][0] << " " << arr[i][1] << " " << arr[i][2] << " " << U(arr[i][0], arr[i][1], arr[i][2]) << endl;
	//}
	in.close();
	out.close();
}

void zeroFirstNValues(std::map<double, double>& myMap, size_t n) {
	auto it = myMap.begin();
	for (size_t i = 0; i < n && it != myMap.end(); ++i, ++it) {
		it->second = 11.22;
	}
}


void generate_cross(vector<double>& gridx, vector<double>& gridy, vector<double>& gridz)
{
	ofstream out("cross.txt");
	int n = gridx.size() * gridy.size() * gridz.size();
	//out << n << endl;
	out << n << endl;
	//double length=5, width=4, height=3;//минимальные размеры горки(ямы?) по x,y,z
	//облдасть горки будет [-length/2;lengrh/2]*[-width/2;width/2]*[-height/2;height/2]
	double R = 100, r = R, height = 2;//радиус ,высота горки
	double h0 = height, rate = 0.99;//начальный шаг, коэффициент разрядки
	set<double, greater<double>>ring;//сформирует массив r-колец
	ring.insert(r);
	int ind = 0;
	for (int i = 0;i < gridx.size();i++)
	{
		for (int j = 0;j < gridy.size();j++)
		{
			double res = sqrt(pow(gridx[i], 2) + pow(gridy[j], 2));
			if (res <= r)//если попали в область горки, добавляем конкретный радиус для кольца
			{
				ring.insert(res);
			}
		}
	}
	for (auto i : ring)  cout << i << " ";
	//vector<pair<double, double>> dict;
	map<double, double> dict;
	double h = h0;
	cout << endl;
	for (auto i : ring)
	{
		//pair<double, double> p = { i,h };
		dict[i] = h;
		h += 0.02;
		//h /= rate;
	}
	//zeroFirstNValues(dict, 20);
	//горка будет полусферой
	cout << "Gorka dict:\n";
	for (auto i : dict)
	{
		cout << i.first << " " << i.second << endl;
	}
	cout << "End Gorka dict.\n";
	int i1 = 0, i2 = 0, i3 = 0;
	for (int i = 0;i < gridz.size();i++)
	{
		for (int j = 0;j < gridy.size();j++)
		{
			for (int k = 0;k < gridx.size();k++)
			{
				if (i == gridz.size() - 1)//если поверхность
				{
					double res = sqrt(pow(gridx[k], 2) + pow(gridy[j], 2));
					if (res <= r)//если попали в область горки
					{
						out << gridx[k] << " " << gridy[j] << " " << gridz[i] + dict[res] << endl;//делаем горку
						//out << gridx[k] << " " << gridy[j] << " " << gridz[i] << endl;//если не хотим горку
					}
					else
					{
						out << gridx[k] << " " << gridy[j] << " " << gridz[i] << endl;
					}
				}
				else
				{
					out << gridx[k] << " " << gridy[j] << " " << gridz[i] << endl;
				}

			}
		}
	}
	out.close();

	ifstream in("cross.txt");
	int sizeG;
	in >> sizeG;
	vector<vector<double>> G(sizeG, vector<double>(3, 0));
	vector<double> t(3);
	for (int u = 0;u < sizeG;u++)
	{
		in >> t[0] >> t[1] >> t[2];
		G[u] = t;
	}
	in.close();
	out.open("local.txt");
	int temp = 1,i,j,k;
	int nx = gridx.size(), ny = gridy.size(), nz = gridz.size();
	out << (nx - 1) * (ny - 1) * (nz - 1) << endl;
	//printf("%d\n",l);//in << l << endl;
	//lx=iter
	for (i = 0;i < nz - 1;i++, temp = 1 + i * nx * ny)
		for (j = 0;j < ny - 1;j++, temp++)
			for (k = 0;k < nx - 1;k++, temp++)
			{
				out << temp << " " << temp + 1 << " " << temp + nx << " " << temp + nx + 1 << " " << temp + nx * ny << " " << temp + nx * ny + 1 << " " << temp + nx * ny + nx << " " << temp + nx * ny + nx + 1 << endl;
				if (G[temp - 1][2] < 0.) out << 0.001 << endl;
				if (0 <= G[temp - 1][2] && G[temp - 1][2] < 500) out << 0.02625 << endl;
				if (500 <= G[temp - 1][2] && G[temp - 1][2] < 800) out << 0.07475 << endl;
				if (800 <= G[temp - 1][2]) out << 0.1 << endl;
				//in << 1 << endl;
			}

	in.close();


}

double mapRange(double value, double a, double b, double c, double d) {
	//return c + (d - c) * (value - a) / (b - a);
	return d - (d - c) * (value - a) / (b - a);
}

//double GetHeight(double x, double y, double center_x = 0.0, double center_y = 0.0) {
//	const double baseheight = 1000.0;
//	const double x_abs = std::abs(x-center_x);
//	const double y_abs = std::abs(y-center_y);
//	const double max_dist = std::max(x_abs, y_abs);
//
//	// Центральная горка
//	if (x_abs <= 20.0 && y_abs <= 20.0) {
//		return baseheight + 12.0;
//	}
//	// Первое кольцо (наклон от 12 до 10)
//	else if (x_abs <= 40.0 && y_abs <= 40.0) {
//		return baseheight + 10.0;
//	}
//	// Второе кольцо (наклон от 10 до 8)
//	else if (x_abs <= 60.0 && y_abs <= 60.0) {
//		return baseheight + 8.0;
//	}
//	// Третье кольцо (наклон от 8 до 6)
//	else if (x_abs <= 80.0 && y_abs <= 80.0) {
//		return baseheight + 6.0;
//	}
//	// Четвёртое кольцо (наклон от 6 до 4)
//	else if (x_abs <= 100.0 && y_abs <= 100.0) {
//		return baseheight + 4.0;
//	}
//	// Пятое кольцо (наклон от 4 до 0)
//	else if (x_abs <= 250.0 && y_abs <= 250.0) {
//		return baseheight + 2.0;
//	}
//	// Шестое кольцо
//	//else if (x_abs <= 450.0 && y_abs <= 450.0) {
//	//	return baseheight + 1.0;
//	//}
//	// Вне горки — оставляем как есть
//	else {
//		return baseheight;
//	}
//}

double GetHeight(
	double x,
	double y,
	double center_x = 0.0,
	double center_y = 0.0,
	double height_hill = 12.0,   // Максимальная высота холма
	double hill_width = 500.0,   // Ширина холма (по X)
	double hill_height = 500.0,  // Высота холма (по Y)
	int steps = 5                // Количество ступенек
) {
	const double baseheight = 1000.0;

	// Смещаем координаты относительно центра
	double dx = x - center_x;
	double dy = y - center_y;

	// Проверяем, что точка внутри прямоугольной области холма
	bool is_inside_hill = (std::abs(dx) <= hill_width / 2.0) && (std::abs(dy) <= hill_height / 2.0);
	if (!is_inside_hill) {
		return baseheight;
	}

	// Нормализуем расстояния от центра до краёв по X и Y
	double nx = std::abs(dx) / (hill_width / 2.0);   // от 0 (центр) до 1 (край)
	double ny = std::abs(dy) / (hill_height / 2.0);

	// Комбинированный коэффициент (чебышёвское расстояние)
	double distance_factor = std::max(nx, ny);

	// Параболическая высота (1 - x²)
	double parabolic_height = height_hill * (1.0 - distance_factor * distance_factor);

	// Дискретизируем высоту (ступеньки) с гарантией достижения пика
	double step_size = height_hill / steps;
	double discretized_height;

	if (distance_factor <= std::numeric_limits<double>::epsilon()) {
		// В центре холма — сразу возвращаем пиковую высоту
		discretized_height = height_hill;
	}
	else {
		// Для остальных точек — округляем вниз с добавлением половины шага для коррекции
		discretized_height = std::floor((parabolic_height + step_size / 2.0) / step_size) * step_size;
	}

	// Гарантируем, что высота не ниже базовой и не выше пиковой
	discretized_height = std::clamp(discretized_height, 0.0, height_hill);

	return baseheight + discretized_height;
}

bool IsInsideHillArea(double x, double y) {
	const double hill_x_min = -200.0, hill_x_max = 200.0;
	const double hill_y_min = -200.0, hill_y_max = 200.0;
	return (x >= hill_x_min) && (x <= hill_x_max) &&
		(y >= hill_y_min) && (y <= hill_y_max);
}

double CalculateHillHeight(double x, double y, double surface_z) {
	const double max_hill_height = 50.0;

	if (!IsInsideHillArea(x, y)) {
		return surface_z;
	}

	double x_center = 0.0, y_center = 0.0; // Центр в (0,0)
	double x_norm = x / 200.0; // Нормировка к [-1,1]
	double y_norm = y / 200.0;

	double r_squared = x_norm * x_norm + y_norm * y_norm;
	if (r_squared > 1.0) {
		return surface_z;
	}

	return surface_z + max_hill_height * (1.0 - r_squared);
}



void relief()
{
	vector<double> gx, gy, gz;
	readgrid("gridx.txt", "gridy.txt", "gridz.txt", gx, gy, gz);
	refineGrid(gx, 1);
	refineGrid(gy, 1);
	refineGrid(gz, 1);
	double r = 10;
	ifstream in;
	
	const double surface_z = gz.back();
	ofstream out("cross.txt");
	ofstream out2("hill.txt");
	int n = gx.size() * gy.size() * gz.size();
	out << n;
	out << std::fixed << setprecision(16) << endl;
	int i1 = 0, i2 = 0, i3 = 0;
	// Генерация сетки
	for (const auto z : gz) {
		for (const auto y : gy) {
			for (const auto x : gx) {
				double out_z = z;

				//if (IsInsideHillArea(x, y))
				//{
				//	
				//	double h = CalculateHillHeight(x, y, surface_z)-1000.0;//1000-998=2
				//	if (z + h >= -h) out_z = z - h;//out << x << " " << y << " " << z - h << " " << endl;
				//
				//}
				//bool is_surface = (z == surface_z);
				//if (is_surface) 
				//{
				//	out2 << x << " " << y << " " << out_z << "\n";
				//	}
				//// Все точки записываем в cross.txt
				//out << x << " " << y << " " << out_z << "\n";


				bool is_surface = (z == surface_z);
				if (is_surface) {
					//out_z = CalculateHillHeight(x, y, surface_z);
					// Записываем в hill.txt только поверхностные точки
					out2 << x << " " << y << " " << out_z << "\n";
				}
				// Все точки записываем в cross.txt
				out << x << " " << y << " " << out_z << "\n";
			}
		}
	}

	std::cout << "Successfully generated:\n";
	out.close();

	in.open("cross.txt");
	int sizeG;
	in >> sizeG;
	vector<vector<double>> G(sizeG, vector<double>(3, 0));
	vector<double> t(3);
	for (int u = 0;u < sizeG;u++)
	{
		in >> t[0] >> t[1] >> t[2];
		G[u] = t;
	}
	in.close();
	out.open("local.txt");
	int temp = 1, i, j, k;
	int nx = gx.size(), ny = gy.size(), nz = gz.size();
	out << (nx - 1) * (ny - 1) * (nz - 1) << endl;
	//printf("%d\n",l);//in << l << endl;
	//lx=iter
	for (i = 0;i < nz - 1;i++, temp = 1 + i * nx * ny)
		for (j = 0;j < ny - 1;j++, temp++)
			for (k = 0;k < nx - 1;k++, temp++)
			{
				out << temp << " " << temp + 1 << " " << temp + nx << " " << temp + nx + 1 << " " << temp + nx * ny << " " << temp + nx * ny + 1 << " " << temp + nx * ny + nx << " " << temp + nx * ny + nx + 1 << endl;
				//if (G[temp - 1][2] < 0.) out << 0.001 << endl;
				//if (0 <= G[temp - 1][2] && G[temp - 1][2] < 500) out << 0.02625 << endl;
				//if (500 <= G[temp - 1][2] && G[temp - 1][2] < 800) out << 0.07475 << endl;
				//if (800 <= G[temp - 1][2]) out << 0.1 << endl;
				out << 1 << endl;
			}

	in.close();

	
	ofstream in2("kraev_1.txt");
	int a, b, c, d;

	in2 << std::fixed << setprecision(16) << endl;

	//нижн¤¤ грань 1
	temp = 1;
	for (i = 0;i < nx * ny;i++, temp++)
	{
		in2 << temp << " " << ug(G[temp - 1]) << endl;
	}
	in2 << endl;
	//лева¤ грань 2
	temp = 1;
	for (i = 0;i < ny * nz;i++)
	{
		in2 << temp << " " << ug(G[temp - 1]) << endl;
		temp += nx;
	}
	in2 << endl;
	//передн¤¤ грань 3
	int g = 1;
	temp = 1;
	for (i = 0;i < nx * nz;i++)
	{
		in2 << temp << " " << ug(G[temp - 1]) << endl;
		if (g == nx)
		{
			temp += nx * ny - nx + 1;
			g = 1;
		}
		else
		{
			temp++;
			g++;
		}
	}
	in2 << endl;
	//верхн¤¤ 4
	temp = 1 + nx * ny * (nz - 1);
	for (i = 0;i < nx * ny;i++)
	{
	    in2 << temp << " " << ug(G[temp - 1]) << endl;
	    temp++;
	}
	in2 << endl;
	//права¤ грань
	temp = nx;
	for (i = 0;i < ny * nz;i++)
	{
		in2 << temp << " " << ug(G[temp - 1]) << endl;
		temp += nx;
	}
	in2 << endl;
	//задн¤¤ грань
	temp = nx * ny - nx + 1;
	g = 1;
	for (i = 0;i < nx * nz;i++)
	{
		in2 << temp << " " << ug(G[temp - 1]) << endl;
		if (g == nx) {
			temp += nx * ny - nx + 1;
			g = 1;
		}
		else
		{
			temp++;
			g++;
		}
	}
	in2.close();


}

int main()
{
	ifstream in;
	vector<double>gridx, gridy, gridz;
	//readgrid("gridx.txt", "gridy.txt", "gridz.txt", gridx, gridy, gridz);
	//generate_cross(gridx, gridy, gridz);
	relief();
	//gen_file("cross.txt","local.txt",gridx, gridy, gridz, true);//основная сетка
	//gen_file("cross_par.txt", "local_par.txt", gridx, gridy, gridz, false);//сетка с переинтерполяцией


	ofstream out("receiver.txt");
	out << std::fixed << setprecision(16) << endl;
	GLOBAL_MATRIX my;
	my.kraev_1();
	my.MSG();
	my.print_result();
	//my.Receiver(out, 100 + 400, -290, 1000);
	out.close();
	//my.print_receiver();
	return 0;
}