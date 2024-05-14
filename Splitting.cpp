#include <iostream>
#include <cmath>
#include <vector>

#include "Task.hpp"

namespace Progonka
{
    double gx(double x)
    {
        return 1 + x;
    }

    double px(double x)
    {
        return 1 + (x * x);
    }

    double ux(double x)// функция u(x)
    {
        return (x * x) * ((1 - x) * (1 - x));
    }

    double fx(double x)// функция f(x), вычисляется по формуле с двумя производными
    {
        return  -2 + 12 * x - 17 * (x * x) + 23 * (x * x * x) - 21 * std::pow(x, 4)  + std::pow(x, 5);
    }

    int n = 10; // кол-во делений
    double h = 1.0 / n;// шаг, используется для итерации

    double* ai = new double[n + 1];
    double* qi = new double[n + 1];
    double* fi = new double[n + 1];
    double* ui = new double[n + 1];
    double* alpha = new double[n + 1];
    double* betta = new double[n + 1];
    double* a = new double[n + 1];
    double* b = new double[n + 1];
    double* c = new double[n + 1];
    double* ypi = new double[n + 1];
    double* yi = new double[n + 1];
    double* y_prev = new double[n + 1];

    void set()
    {
        ypi[0] = 0;
        ypi[n] = 0;

        for (int i = 0; i <= n; i++) // заполнение массивов, вместо u(x) будет u(ih) и т.д.
        {
            ai[i] = px(i * h);
            qi[i] = gx(i * h);
            fi[i] = fx(i * h) * h * h;
            ui[i] = ux(i * h);
        }

        for (int i = 0; i < n; i++)// вычисление коэффициентов A,B,C - нужны для вычисления alpha,betta по формулам метода прогонки
        {
            a[i] = -ai[i];
            b[i] = -(ai[i] + ai[i + 1] + h * h * qi[i]);
            c[i] = -ai[i + 1];
        }

        alpha[0] = c[0] / b[0];
        betta[0] = fi[0] / b[0];

        for (int i = 0; i < n; i++)// вычисление коэффициентов alpha,betta
        {
            alpha[i + 1] = c[i] / (b[i] - (a[i] * alpha[i]));
            betta[i + 1] = ((a[i] * betta[i]) - fi[i]) / (b[i] - (a[i] * alpha[i]));
        }

        for (int i = 0; i <= n; ++i) {
            yi[i] = 0; 
            y_prev[i] = yi[i];
        }
    }

    void main1()
    {
        for (int i = n - 1; i >= 0; i--)// вычисление результата(yi - вектор решений)
        {
            ypi[i] = alpha[i + 1] * ypi[i + 1] + betta[i + 1];
        }

        ypi[n] =  betta[n + 1];

        //вывод таблицы
        std::cout << "i*h" << "     y[i]" << "            u(ih)" << "         |y[i]-U(ih)|" << std::endl;
        for (int i = 0; i <= n; i++)
        {
            std::cout << i * h << "\t"; //узлы
            std::cout << ypi[i] << "\t";//соответствующие значения из вектора
            std::cout << ui[i] << "\t";//значение функции u(x) в узлах
            std::cout << std::abs(ypi[i] - ui[i]) << std::endl;// погрешность решения
        }

        delete[] ai;
        delete[] qi;
        delete[] fi;
        delete[] ui;
        delete[] alpha;
        delete[] betta;
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] ypi;
        delete[] yi;
        delete[] y_prev;
    }

    

    void main2()
    {  
        //начальные значения
        yi[0] = 0;
        yi[n] = 0;

        for (int i = n - 1; i >= 0; i--)// вычисление результата(yi - вектор решений)
        {
            ypi[i] = alpha[i + 1] * ypi[i + 1] + betta[i + 1];
        }

        ypi[n] =  betta[n + 1];

        double eps = 1e-6; // точность
        double error = eps + 1; // начальное значение ошибки
        int iterations = 0; // количество итераций

        // Итерационный процесс
        while (error > eps && iterations < 1000)
        {
            iterations++;
            error = 0;
            yi[0] = (ai[1]*y_prev[1] + fi[0])/-b[0];
            for (int i = 1; i < n; i++)
            {
                yi[i] = (ai[i]*yi[i-1] + ai[i+1]*y_prev[i+1] + fi[i])/-b[i];
            }
            //yi[n] = (ai[n]*yi[n-1] + fi[n])/-b[n];
            for (int i = 0; i <= n; i++)
                error = std::max(error, a[i]);

            std::copy(yi, yi + n + 1, y_prev); // Обновляем вектор приближения
             
        }

        //вывод таблицы
        std::cout << "i*h" << "     yp[i]" << "\t\ty[i]" << "\t\tпогр" << std::endl;
        for (int i = 0; i <= n; i++)
        {
            std::cout << i * h << "\t"; //узлы
            std::cout << ypi[i] << "\t";
            std::cout << yi[i] << "\t";//соответствующие значения из вектора
            std::cout << std::abs(ypi[i] - yi[i]) << "\t" << std::endl;
        }

        std::cout << "Iterations: " << iterations << std::endl;

        delete[] ai;
        delete[] qi;
        delete[] fi;
        delete[] ui;
        delete[] alpha;
        delete[] betta;
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] yi;
        delete[] y_prev;
    }

    void main3()
    {
        //начальные значения
        yi[0] = 0;
        yi[n] = 0;

        double eps = 1e-6; // точность
        //double error = eps + 1; // начальное значение ошибки
        int iterations = 0; // количество итераций

        double error = std::numeric_limits<double>::max(); // начальное значение ошибки

        for(double w = 1; w > eps; w -= h)
        {
            iterations = 0; // обнуляем количество итераций перед каждой итерацией с новым w
            error = std::numeric_limits<double>::max(); // обнуляем ошибку перед каждой итерацией с новым w

            while (error > eps && iterations < 1000)
            {
                iterations++;
                error = 0;
                yi[0] = (1-w)*(y_prev[1]) + w*((ai[1]*y_prev[1] + fi[0])/-b[0]);
                for (int i = 1; i < n; i++)
                {
                    yi[i] = (1-w)*(y_prev[i+1]) + w*((ai[i]*yi[i-1] + ai[i+1]*y_prev[i+1] + fi[i])/-b[i]);
                }

                for (int i = 1; i < n; i++)
                    error = std::max(error, std::abs(yi[i] - y_prev[i]));

                std::copy(yi, yi + n + 1, y_prev); // Обновляем вектор приближения
            }

            std::cout << w << "\t"; //узлы
            std::cout << iterations << "\t" << std::endl;//соответствующие значения из вектора
        }


        delete[] ai;
        delete[] qi;
        delete[] fi;
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] yi;
        delete[] y_prev;
    }

}

int main()
{
    Progonka::set();
    Progonka::main3();


  return 0;

}
