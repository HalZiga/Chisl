#include <iostream>
#include <cmath>

#define M_PI 3.14159265358979323846

static const double a = 0.0;
static const double b = 2.0;
static const double eps = 0.000001;

double MyErf(double x) {
    double sum, n, a, q;
    n = 0;
    a = (2 * x)/(sqrt(M_PI));
    q = 0;
    sum = a;
    while (fabs(a) > 1e-6)
    {
        q = (-1) * ((((2 * n) + 1)) / ((2 * (n * n)) + (5 * n) + 3)) * (x * x);
        a *= q;
        sum += a;
        n++;
    }
    return sum;
}

double Podint(double t) {
    return std::exp(-1 * t * t) * 2 / std::sqrt(M_PI);
}

double RightRectangles(double x, int& n) {
    double s1 = 0, s = 0;
    n = 1;
    double h;
    do {
        s = s1;
        n *= 2;
        h = (x - a) / n;
        s1 = 0;
        for (double tempX = 0; tempX < x; tempX += h) {
            s1 += h * (Podint(tempX + h));
        }
    } while (std::abs(s1 - s) > eps);
    return (s1);
}

double CentralRectangles(double x, int& n) {
    double s1 = 0, s = 0;
    n = 1;
    double h;
    do {
        s = s1;
        n *= 2;
        h = (x - a) / n;
        s1 = 0;
        for (double tempX = 0; tempX < x; tempX += h) {
            s1 += h * (Podint(tempX + h / 2));
        }
    } while (std::abs(s1 - s) > eps);
    return s1;
}

double Trapeciya(double x, int& n) {
    double s1 = 0, s = 0;
    n = 1;
    double h;
    do {
        s = s1;
        n *= 2;
        h = (x - a) / n;
        s1 = 0;
        for (double tempX = 0; tempX < x; tempX += h) {
            s1 += h * ((Podint(tempX) + Podint(tempX + h)) / 2);
        }
    } while (std::abs(s1 - s) > eps);
    return s1;
}

double Simpson(double x, int& n) {
    double s1 = 0, s = 0;
    n = 1;
    double h;
    do {
        s = s1;
        n *= 2;
        h = (x - a) / n;
        s1 = 0;
        for (double tempX = 0; tempX < x; tempX += h) {
            s1 += h / 6 * (Podint(tempX) + 4 * Podint(tempX + h / 2) + Podint(tempX + h));
        }
    } while (std::abs(s1 - s) > eps);
    return s1;
}

double Gauss(double x, int& n) {
    double s1 = 0, s;
    n = 1;
    double h;
    do {
        s = s1;
        n *= 2;
        h = (x - a) / n;
        s1 = 0;
        for (double tempX = 0; tempX < x; tempX += h) {
            s1 += h / 2 * (Podint(tempX + (h / 2) * (1 - 1 / std::sqrt(3))) +
                           Podint(tempX + (h / 2) * (1 + 1 / std::sqrt(3))));
        }
    } while (std::abs(s1 - s) > eps);
    return s1;
}

int main() {
    double h = 0.2;
    int n = 1;

    std::cout << "Метод правых прямоугольников\n";
    std::cout << "xi\tJ_0(x)\t\tJ_N(x)\t\t|J_0(x)-J_N(x)|\tN\n";
    for (int i = 0; i * h <= b; i++) {
        std::cout << std::round(i * h * 10) / 10.0 << "\t" << MyErf(i * h) << "\t" << RightRectangles(i * h, n) << "\t" << std::abs(MyErf(i * h) - RightRectangles(i * h, n)) << "\t" << n << "\n";
    }

    std::cout << "\nМетод центральных прямоугольников\n";
    std::cout << "xi\tJ_0(x)\t\tJ_N(x)\t\t|J_0(x)-J_N(x)|\tN\n";
    for (int i = 0; i * h <= b; i++) {
        std::cout << std::round(i * h * 10) / 10.0 << "\t" << MyErf(i * h) << "\t" << CentralRectangles(i * h, n) << "\t" << std::abs(MyErf(i * h) - CentralRectangles(i * h, n)) << "\t" << n << "\n";
    }

    std::cout << "\nМетод трапеций\n";
    std::cout << "xi\tJ_0(x)\t\tJ_N(x)\t\t|J_0(x)-J_N(x)|\tN\n";
    for (int i = 0; i * h <= b; i++) {
        std::cout << std::round(i * h * 10) / 10.0 << "\t" << MyErf(i * h) << "\t" << Trapeciya(i * h, n) << "\t" << std::abs(MyErf(i * h) - Trapeciya(i * h, n)) << "\t" << n << "\n";
    }

    std::cout << "\nМетод Симпсона\n";
    std::cout << "xi\tJ_0(x)\t\tJ_N(x)\t\t|J_0(x)-J_N(x)|\tN\n";
    for (int i = 0; i * h <= b; i++) {
        std::cout << std::round(i * h * 10) / 10.0 << "\t" << MyErf(i * h) << "\t" << Simpson(i * h, n) << "\t" << std::abs(MyErf(i * h) - Simpson(i * h, n)) << "\t" << n << "\n";
    }

    std::cout << "\nМетод Гаусса\n";
    std::cout << "xi\tJ_0(x)\t\tJ_N(x)\t\t|J_0(x)-J_N(x)|\tN\n";
    for (int i = 0; i * h <= b; i++) {
        std::cout << std::round(i * h * 10) / 10.0 << "\t" << MyErf(i * h) << "\t" << Gauss(i * h, n) << "\t" << std::abs(MyErf(i * h) - Gauss(i * h, n)) << "\t" << n << "\n";
    }

    return 0;
}