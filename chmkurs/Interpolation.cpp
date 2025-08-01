#include "myclass.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ranges>
#include <iomanip>
#include <chrono>  
#include <omp.h>
#include <numeric>
using namespace std;

//double ug(vector<double> a)
//{
//    return 0;
//    //return a[0]+a[1]+a[2];
//}


// Константы квадратуры Гаусса
//constexpr array<double, 3> tau = { 8. / 9., 5. / 9., 5. / 9. };
//constexpr array<double, 3> t = { 0., 0.77459666924148337, -0.77459666924148337 };

double basic_func(int i, double x, double y, double z)//для единичного куба
{
    double X[2] = { 1 - x,x };
    double Y[2] = { 1 - y,y };
    double Z[2] = { 1 - z,z };

    switch (i)
    {
    case 0: return X[0] * Y[0] * Z[0];
    case 1: return X[1] * Y[0] * Z[0];
    case 2: return X[0] * Y[1] * Z[0];
    case 3: return X[1] * Y[1] * Z[0];
    case 4: return X[0] * Y[0] * Z[1];
    case 5: return X[1] * Y[0] * Z[1];
    case 6: return X[0] * Y[1] * Z[1];
    case 7: return X[1] * Y[1] * Z[1];
    }
}


double basic_func(double x0, double x1, double y0, double y1, double z0, double z1, int i, double x, double y, double z)
{
    double xi = (x - x0) / (x1 - x0);
    double eta = (y - y0) / (y1 - y0);
    double zeta = (z - z0) / (z1 - z0);

    double X[2] = { 1 - xi, xi };
    double Y[2] = { 1 - eta, eta };
    double Z[2] = { 1 - zeta, zeta };

    int ix = i % 2;
    int iy = (i / 2) % 2;
    int iz = i / 4;
    return X[ix] * Y[iy] * Z[iz];
}

struct Point
{
    double x, y, z;
    Point(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    Point()
    {
        x = 0.;
        y = 0.;
        z = 0.;
    }
};

//void print(vector<double> a)
//{
//    for (int i = 0;i < a.size();i++)
//        cout << a[i] << " ";
//    cout << endl;
//}


void findEL(vector<double>x, vector<double>y, vector<double> z, double h)
{
    int  i, j, k, flag = 3;
    vector<double> t1, t2, t3;
    for (i = 1;i < x.size() - 1 && flag;i++)
    {
        if (x[i] - x[i - 1] == h && x[i + 1] - x[i] == h) {
            t1.push_back(x[i]);
            flag--;
        }
    }
    flag = 3;
    for (i = 1;i < y.size() - 1 && flag;i++)
    {
        if (y[i] - y[i - 1] == h && y[i + 1] - y[i] == h) {
            t2.push_back(y[i]);
            flag--;
        }
    }
    flag = 3;
    for (i = 1;i < z.size() - 1 && flag;i++)
    {
        if (z[i] - z[i - 1] == h && z[i + 1] - z[i] == h) {
            t3.push_back(z[i]);
            flag--;
        }
    }
    t3.push_back(z[z.size() - 1]);
    ofstream in("../../chmkurs/chmkurs/el.txt");
    in << std::fixed << setprecision(16) << endl;
    in << t1[0] << " " << t2[0] << " " << t3[1] << endl;//Общий
    in << t1[1] << " " << t2[0] << " " << t3[1] << endl;//правый
    in << t1[0] << " " << t2[0] << " " << t3[0] << endl;//нижний
    // cout  << t1 << " " << t2 << " " << t3;
    in.close();
}



// Рекурсивная функция построения KD-Tree
vector<size_t> build_kdtree(vector<ElementBounds>& elements) {
    vector<size_t> indices(elements.size());
    iota(indices.begin(), indices.end(), 0);

    // Быстрая сортировка по осям X/Y/Z без рекурсии
    auto sort_axis = [&](int axis) {
        sort(indices.begin(), indices.end(),
            [&elements, axis](size_t a, size_t b) {
                return axis == 0 ? (elements[a].x_min < elements[b].x_min) :
                    axis == 1 ? (elements[a].y_min < elements[b].y_min) :
                    (elements[a].z_min < elements[b].z_min);
            });
        };

    // Чередуем оси для балансировки
    sort_axis(0);
    sort_axis(1);
    sort_axis(2);

    return indices;
}


//для отчета первые тесты были с numvel=numgel=10
//во втором исследовании numgel=10,numvel=7

GLOBAL_MATRIX::GLOBAL_MATRIX(vector<cross>& Cross,vector<local>& Matr,vector<double>& q)//входные параметры - это из шестигранника
{

    
    
    
    ifstream in;

    read_cross("cross_par.txt");//прямоугольный
    read_local("local_par.txt");// прямоугольный

    vector<ElementBounds> elements;
    for (size_t i = 0; i < matr.size(); i++) {
        // Вычисляем границы по ВСЕМ 8 вершинам (важно для регулярных сеток!)
        double x_min = set[matr[i].mas[0] - 1].mas[0], x_max = x_min;
        double y_min = set[matr[i].mas[0] - 1].mas[1], y_max = y_min;
        double z_min = set[matr[i].mas[0] - 1].mas[2], z_max = z_min;

        for (int j = 1; j < 8; j++) {
            int v_idx = matr[i].mas[j] - 1;
            x_min = min(x_min, set[v_idx].mas[0]);
            x_max = max(x_max, set[v_idx].mas[0]);
            y_min = min(y_min, set[v_idx].mas[1]);
            y_max = max(y_max, set[v_idx].mas[1]);
            z_min = min(z_min, set[v_idx].mas[2]);
            z_max = max(z_max, set[v_idx].mas[2]);
        }
        elements.push_back({ x_min, x_max, y_min, y_max, z_min, z_max, (int)i });
    }
    
   
    vector<size_t> kdtree_order = build_kdtree(elements);

    auto kdtree_search = [&](double x, double y, double z) -> int {
        for (size_t idx : kdtree_order) {
            const auto& e = elements[idx];
            // Ранний выход при первой же найденной точке
            if (x >= e.x_min && x <= e.x_max &&
                y >= e.y_min && y <= e.y_max &&
                z >= e.z_min && z <= e.z_max) {
                return e.original_index;
            }
        }
        return -1;
        };

    d.resize(n);
    f.resize(n);
    x.resize(n);
    ig.resize(n + 1);

    formier_profil();//сформируем профиль для регулярной сетки(set и matr прямоугольные)
    vector<vector<PointData>> fe(matr.size());//сделаем массив,хранящий по индексу параллелепипеда точки Гаусса с их значениями
    bool flag = false;
   
    //заполним этот массив при помощи прохода по шестигранникам, каждую точку будем по её локальным координатам находить глобальные координаты(долго)
    cout << "Filling parallelepipeds with Gaussian points\n";
    double t[3] = {0., 0.77459666924148337, -0.77459666924148337};
    double xg=0., yg=0., zg=0.,res=0.;

    vector<vector<double>> basis(8, vector<double>(27));
    for (int gi = 0; gi < 3; ++gi) {
        double theta = 0.5 + t[gi] * 0.5;
        for (int gj = 0; gj < 3; ++gj) {
            double eta = 0.5 + t[gj] * 0.5;
            for (int gk = 0; gk < 3; ++gk) {
                double xi = 0.5 + t[gk] * 0.5;
                int idx = gi * 9 + gj * 3 + gk;
                for (int i = 0; i < 8; ++i) {
                    basis[i][idx] = basic_func(i, xi, eta, theta);
                }
            }
        }
    }
    

    for (int el = 0;el < Matr.size();el++)
    {
        //проходимся по всем 27 точкам Гаусса
        for (int gi = 0; gi < 3; ++gi) {
            for (int gj = 0; gj < 3; ++gj) {
                for (int gk = 0; gk < 3; ++gk) {
                    int idx = gi * 9 + gj * 3 + gk;
                    xg = 0;
                    yg = 0.;
                    zg = 0.;
                    res = 0.;
                    //определяем глобальные координаты гауссовой точки и её значение
                    for (int i = 0;i < 8;i++)
                    {
                        xg += Cross[Matr[el].mas[i]-1].mas[0]*basis[i][idx];
                        yg += Cross[Matr[el].mas[i]-1].mas[1]* basis[i][idx];
                        zg += Cross[Matr[el].mas[i]-1].mas[2]* basis[i][idx];
                        res += q[Matr[el].mas[i] - 1] * basis[i][idx];
                    }
                    int i;
                    flag = false;
                    
                    int found_idx = kdtree_search(xg, yg, zg);
                        //находим нужный элемент в параллелепипедной сетке
                        //for (i = 0; i < matr.size() && !flag; i++) {
                        //    int v1 = matr[i].mas[0];
                        //    int v2 = matr[i].mas[1];
                        //    int v3 = matr[i].mas[2];
                        //    int v5 = matr[i].mas[4];
                        //    //проверка, чтобы в i-ом КЭ наша точка лежала внутри него
                        //    if (xg >= set[v1 - 1].mas[0] && xg <= set[v2 - 1].mas[0] &&
                        //        yg >= set[v1 - 1].mas[1] && yg <= set[v3 - 1].mas[1] &&
                        //        zg >= set[v1 - 1].mas[2] && zg <= set[v5 - 1].mas[2])
                        //        flag = true;
                        //}
                        //cout << found_idx << endl;
                        //заносим в i параллелепипед гауссову точку
                        if (found_idx!=-1)
                            fe[found_idx].push_back(PointData(xg, yg, zg, res));
                    
                }
            }
        }
        
    }
    cout << "Filling done!\n";
    cout << "Formier Matrix the least squares method\n";
    //после прохода по шестигранникам будет проход по параллелепипедам(массив уже готов потому что)
    for (int el = 0;el < matr.size();el++)
    {
    //в этом проходе будем генерировать для каждого КЭ локальную матрицу и вектор по формулам из МНК(долго)
        vector<vector<double>> A(8, vector<double>(8, 0.0));
        vector<double> b(8, 0);
        double x0 = set[matr[el].mas[0] - 1].mas[0];
        double y0 = set[matr[el].mas[0] - 1].mas[1];
        double z0 = set[matr[el].mas[0] - 1].mas[2];

        double x1 = set[matr[el].mas[1] - 1].mas[0];
        double y1 = set[matr[el].mas[2] - 1].mas[1];
        double z1 = set[matr[el].mas[4] - 1].mas[2];
        for (int i = 0;i < 8;i++)
        {
            for (int j = 0;j < 8;j++)
            {
                for (int k = 0;k < fe[el].size();k++)
                {
                    
                    A[i][j] += basic_func(x0,x1,y0,y1,z0,z1,i,fe[el][k].x,fe[el][k].y,fe[el][k].z)
                        * basic_func(x0, x1, y0, y1, z0, z1, j, fe[el][k].x, fe[el][k].y, fe[el][k].z);
                }
            }
            for (int k = 0;k < fe[el].size();k++)
            {
                b[i] += basic_func(x0,x1,y0,y1,z0,z1,i, fe[el][k].x, fe[el][k].y, fe[el][k].z)
                    * fe[el][k].value;
            }
        }
        //затем локальную матрицу вносить в глобальную так же, как при сборке матрицы жесткости
        
        for (int i = 0; i < 8; i++) {
            int global_i = matr[el].mas[i] - 1;
            d[global_i] += A[i][i];
            //arr[global_i][global_i] += A[i][i];
            f[global_i] += b[i];
            for (int j = 0; j < i; j++) {
                int global_j = matr[el].mas[j] - 1;
                double value = A[i][j];

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
    cout << "Matrix formier!\n";
    cout << "MSG...\n";
    MSG();
    cout << "MSG complete!\n";
    for (int i = 0;i < 10;i++)    cout << x[i] << " ";
    cout << endl;
    //разложение в интересующих точках по весам и базисным параллелепипедам, получаем решение(быстро)
    ofstream out("interp.txt");
    vector<double> a;
    for (int i = 0;i < 69;i++)
    {
        out << set[i].mas[0] << " " << U(set[i].mas[0], 0, 1000)<<endl;
    }
    cout << "Value: " << U(1.5,1,3) << endl;

    
}

