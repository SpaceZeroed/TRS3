// TRS3.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
#define PI acos(-1.)
using namespace std;
namespace var9
{
    double True_U(double x, double y)
    {
        return sin( PI * x * x) * sin(PI * y / 2);
    }
    double F(double x, double y)
    {
        return 2 * PI * cos(PI * x * x) * sin(PI * y / 2)
            - 4 * PI * PI * x * x * sin(PI * x * x) * sin(PI * y / 2)
            - sin(PI * x * x) * PI * PI * sin(PI * y / 2) / 4;
    }
    double Lx = 1, Ly = 2;
    double phi_1(double y)
    {
        return 0;
    }
    double phi_2(double y)
    {
        return sin(PI * Lx * Lx) * sin(PI * y / 2);
    }
    double phi_3(double x)
    {
        return 0;
    }
    double phi_4(double x)
    {
        return sin(PI * x * x) * sin(PI * Ly / 2);
    }
};
void PrintMatrix(vector<vector<double>> Matrix)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            cout << setw(7) << Matrix[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
vector<double> Matrix_Vector_Multiplication(int n, int m , double alpha, double betta, double gamma, vector<double> f)
{
    double new_betta = betta / alpha, new_gamma = gamma / alpha;
    vector<double> temp(f.size());
    // задал решение для первой строчки u 
    temp[0] = f[1] * new_betta + f[n - 1] * new_gamma;
    for (int i = 1; i < n - 2; i++)
    {
        temp[i] = (f[i + 1] + f[i - 1]) * new_betta + f[i + (n - 1)] * new_gamma;
    }
    temp[n - 2] = f[n - 3] * new_betta + f[n - 2 + n - 1] * new_gamma;
    // задаю значения для середины 
    for (int i = 1; i <= m - 3; i++)
    {
        temp[i * (n - 1)] = f[i * (n - 1) + 1] * new_betta + ( f[i * (n - 1) + (n - 1)] + f[i * (n - 1) - (n - 1)] ) * new_gamma;
        for (int j = 1; j < n - 2; j++)
        {
            temp[i * (n - 1) + j] = ( f[i * (n - 1) + j + 1] + f[i * (n - 1) + j - 1]) * new_betta + ( f[i * (n - 1) + j + (n - 1)] + f[i * (n - 1) + j - (n - 1)] ) * new_gamma;
        }
        temp[i * (n - 1) + (n - 2)] = f[i * (n - 1) + (n - 3) ] * new_betta + (f[i * (n - 1) + (n - 2) + (n - 1)] + f[i * (n - 1) + (n - 2) - (n - 1)]) * new_gamma;
    }
    // задаю решение для последней строчки 
    temp[ (m - 2) * (n - 1)] = f[1] * new_betta + f[(m - 2) * (n - 1) - (n - 1)] * new_gamma;
    for (int i = (m - 2) * (n - 1) + 1; i < (m - 1) * (n - 1) - 1; i++)
    {
        temp[i] = (f[ i + 1] + f[i - 1]) * new_betta + f[i + (n - 1)] * new_gamma;
    }
    temp[(m - 1) * (n - 1) - 1] = f[(m - 1) * (n - 1) - 2] * new_betta + f[(m - 1) * (n - 1) - 1 - (n - 1)] * new_gamma;
    return temp;
}
using namespace var9;
vector<double> Ex1(int n, int m)
{
    double hx = Lx / n, hy = Ly / m;
    vector<double> X(n + 1);
    vector<double> Y(m + 1);
    for (int i = 0; i <= n; i++)
        X[i] = i * hx;
    for (int i = 0; i <= m; i++)
        Y[i] = i * hy;
    vector<double> f((n - 1) * (m - 1));
    for (int i = 0; i <= m - 2; i++) // у рассматриваем с 1 по m - 1 элемент, а записываем как 0 .. m - 2
    {
        for (int j = 0; j <= n - 2; j++) // х, по аналогии с y с 1 по n - 1 -> 0 .. n - 2
        {
            if (i == 0)
            {
                if (j == 0)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_3(X[j]) / hy / hy - phi_1(Y[i]) / hx / hx ;
                }
                else if( j != n - 2)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_3(X[j]) / hy / hy;
                }
                else
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_3(X[j]) / hy / hy - phi_2(Y[i]) / hx / hx;
                }
            }
            else if (i != 0 && i != m - 2)
            {
                if (j == 0)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_1(Y[i]) / hx / hx;
                }
                else if (j != n - 2)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]);
                }
                else
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_2(Y[i]) / hx / hx;
                }
            }
            else
            {
                if (j == 0)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_4(X[j]) / hy / hy - phi_1(Y[i]) / hx / hx;
                }
                else if (j != n - 2)
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_4(X[j]) / hy / hy;
                }
                else
                {
                    f[i * (n - 1) + j] = F(X[j], Y[i]) - phi_4(X[j]) / hy / hy  - phi_2(Y[i]) / hx / hx;
                }
            }
        }
    }
    double alpha = -2 * (1. / hx / hx + 1. / hy / hy), betta = 1./ hx / hx, gamma = 1./ hy / hy;

    for (int i = 0; i < (n - 1) * (m - 1); i++)
        f[i] = f[i] / alpha;

    vector<double> temp_u = f;
    double raz = 1;
    do
    {
        
    }while (raz > 1e-6)
}
int main()
{
    std::cout << "Hello World!\n";
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
