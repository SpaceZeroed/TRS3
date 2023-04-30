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
        return 2. * PI * cos(PI * x * x) * sin(PI * y / 2.)
            - 4. * PI * PI * x * x * sin(PI * x * x) * sin(PI * y / 2.)
            - sin(PI * x * x) * PI * PI * sin(PI * y / 2) / 4.;
    }
    double Lx = 1, Ly = 2;
    double phi_1(double y)
    {
        return 0;
    }
    double phi_2(double y)
    {
        return 0;
    }
    double phi_3(double x)
    {
        return 0;
    }
    double phi_4(double x)
    {
        return 0;
    }
};
using namespace var9;
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
void PrintVector(vector<double> V)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < V.size(); i++)
    {
        cout << setw(7) << V[i] << "  " << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
vector<double> Make_F_vector(int n, int m)
{
    double hx = Lx / n, hy = Ly / m ;
    vector<double> X(n + 1);
    vector<double> Y(m + 1);
    for (int i = 0; i <= n; i++)
        X[i] = i * hx;
    for (int i = 0; i <= m; i++)
        Y[i] = i * hy;
    vector<double> f((n - 1) * (m - 1));
    //cout << (m - 1) << endl;
    for (int i = 0; i < m - 1; i++) // у рассматриваем с 1 по m - 1 элемент, а записываем как 0 .. m - 2
    {
        for (int j = 0; j < n - 1; j++) // х, по аналогии с y с 1 по n - 1 -> 0 .. n - 2
        {
            if (i == 0)
            {
                if (j == 0)
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_3(X[j + 1 ]) / hy / hy - phi_1(Y[i + 1]) / hx / hx;

                else if (j != n - 2)
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_3(X[j + 1]) / hy / hy;

                else
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_3(X[j + 1]) / hy / hy - phi_2(Y[i + 1]) / hx / hx;

            }
            else if (i != 0 && i != m - 2)
            {
                if (j == 0)
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_1(Y[i + 1]) / hx / hx;

                else if (j != n - 2)
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]);

                else
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_2(Y[i + 1]) / hx / hx;

            }
            else
            {
                if (j == 0)
                {
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_4(X[j + 1]) / hy / hy - phi_1(Y[i + 1]) / hx / hx;
                }
                else if (j != n - 2)
                {
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_4(X[j + 1]) / hy / hy;
                }
                else
                {
                    f[i * (n - 1) + j] = F(X[j + 1], Y[i + 1]) - phi_4(X[j + 1]) / hy / hy - phi_2(Y[i + 1]) / hx / hx;
                }
            }
        }
    }

    //cout << " f size is  " << f.size() << endl;
    return f;
}

vector<double> Matrix_Vector_Multiplication(int n, int m , double alpha, double betta, double gamma, vector<double> f)
{
    double new_betta =  -betta / alpha, new_gamma = -gamma / alpha;
    vector<double> temp(f.size());
    //cout << f.size() << endl;
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
    for (int i = (m - 2) * (n - 1) + 1; i <= (m - 1) * (n - 1) - 2; i++)
    {
        temp[i] = (f[ i + 1] + f[i - 1]) * new_betta + f[i - (n - 1)] * new_gamma;
    }
    temp[(m - 1) * (n - 1) - 1] = f[(m - 1) * (n - 1) - 2] * new_betta + f[(m - 1) * (n - 1) - 1 - (n - 1)] * new_gamma;
    return temp;
}
double NormaRazn(vector<double> a1, vector<double> a2)
{
    double raz = 0;
    for (int i = 0; i < a1.size(); i++)
    {
        raz += (a1[i] - a2[i]) * (a1[i] - a2[i]);
    }
    return raz;
}
double MaxRazn(vector<double> a1, vector<double> a2)
{
    double raz = 0;
    for (int i = 0; i < a1.size(); i++)
    {
        if (abs(a1[i] - a2[i]) > raz)
            raz = abs(a1[i] - a2[i]);
    }
    return raz;
}
vector<double> Ex1(int n, int m)  // n - номер послежней точки, точки идут как 0 .. n 
{
    double hx = Lx / n, hy = Ly /  m;
    double alpha = -2. / hx / hx - 2. / hy / hy, betta = 1. / hx / hx, gamma = 1. / hy / hy;

    vector<double> f = Make_F_vector(n, m);
    for (int i = 0; i < (n - 1) * (m - 1); i++)
        f[i] = f[i] / alpha;
    //PrintVector(f);
    vector<double> temp_u = f;
    //cout << "alpha = " << alpha << ", " << "betta = " << betta << ", " << "gamma = " << gamma << endl;

    double raz = 1;
    do
    {
        vector<double> new_u = Matrix_Vector_Multiplication(n, m, alpha, betta, gamma, temp_u);
        for (int i = 0; i < f.size(); i++)
        {
            new_u[i] = f[i] + new_u[i];
        }
        raz = MaxRazn(new_u, temp_u);
        //cout << raz << endl;
        temp_u = new_u;
    } while (raz > 1e-10);
    //cout << "size is " << temp_u.size() << endl;

    return temp_u;
}
//Метод последовательной верхней релаксации(SOR)
vector<double> SOR(int n, int m, double alpha, double betta, double gamma, vector<double> x, vector<double> b, double w)
{
    vector<double>  x1((n - 1) * (m - 1), 0);
    for (int l = 0; l < (n - 1) * (m - 1); l++)
    {
        double sum1 = 0; double sum2 = 0;
        if (l == 0)
            sum1 += 0;
        else if (l <= (n - 2))
            sum1 += betta * x1[l - 1];
        else if (l % (n - 1) == 0)
            sum1 += gamma * x1[l - (n - 1)];
        else
            sum1 += gamma * x1[l - (n - 1)] + betta * x1[l - 1];
        
        if (l == (n - 1) * (m - 1) - 1 )
            sum1 += alpha * x[l];
        else if ( (n - 1) * (m - 2) <= l)
            sum1 += alpha * x[l] + betta * x[l + 1];
        else if (l % ( n - 1) ==  0)
            sum2 += alpha * x[l] + gamma * x[l + (n - 1)];
        else
            sum1 += alpha * x[l]  + betta * x[l + 1] + gamma * x[l + (n - 1)];
        
        x1[l] = x[l] + w * (b[l] - sum1 - sum2) / alpha;
    }
    
    return x1;
}
vector<double> Ex2(int n, int m) // n - номер послеженй точки 
{
    double hx = Lx / n , hy = Ly / m ;
    double alpha = -2 * (1. / hx / hx + 1. / hy / hy), betta = 1. / hx / hx, gamma = 1. / hy / hy;

    vector<double> f = Make_F_vector(n, m);

    vector<double> u1((n - 1) * (m - 1), 0), u2((n - 1) * (m - 1), 0);
    double raz;

    do
    {
        raz = 0.0;
        u2 = SOR(n, m, alpha, betta, gamma, u1, f, 1);
   
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            if (abs(u2[i] - u1[i]) > raz)
            {
                raz = abs(u2[i] - u1[i]);
            }
        }
        for (int i = 0; i < (n - 1) * (m - 1); i++)
        {
            u1[i] = u2[i];
        }
        //PrintVector(u2);
    } while (raz > 1e-10);

    return u2;
}
vector<double> vector_true_U(int n, int m)
{
    double hx = Lx / n, hy = Ly / m;
    vector<double> X(n + 1);
    vector<double> Y(m + 1);
    for (int i = 0; i <= n; i++)
        X[i] = i * hx;
    for (int i = 0; i <= m; i++)
        Y[i] = i * hy;
    vector<double> temp((n - 1) * (m - 1));
    for (int i = 0; i < m - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            temp[i * (n - 1) + j] = True_U(X[j + 1], Y[i + 1]);
        }
    }
    return temp;
}
vector<double> Ex3(int n, int m)
{
    double hx = Lx / n, hy = Ly / m;
    double alpha = -2 * (1. / hx / hx + 1. / hy / hy), betta = 1. / hx / hx, gamma = 1. / hy / hy;

}
int main()
{
    //Ex1(50, 50);
    //Ex2(100, 100);
    //PrintVector(Ex1(100, 100));
    //PrintVector(Ex2(4, 4));
    cout << "max razn = " << MaxRazn(Ex2(50, 50), vector_true_U(50, 50));
    cout << "max razn = " << MaxRazn(Ex2(50, 50), Ex1(50, 50));
    return 0;
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
