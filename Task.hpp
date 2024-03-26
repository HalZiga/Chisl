#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>

#define M_PI 3.14159265358979323846

int beginPoint, endPoint, Amount_of_points;
float step;

std::vector<std::vector<float>> Experiment_Lagrange_Cheb(std::vector<float> data);
std::vector<std::vector<float>> Experiment_Lagrange_Row(std::vector<float> data);
std::vector<std::vector<float>> Pogreshnost_Lagrange_Row(std::vector<float> data);
std::vector<std::vector<float>> Pogreshnost_Lagrange_Cheb(std::vector<float> data);
float Get_Max_Pogreshnost(std::vector<std::vector<float>> Pogreshnost);
std::vector<std::vector<float>> Get_Pogreshnost(std::vector<std::vector<float>> TabulatePolinom, std::vector<std::vector<float>> InterpolatePolinom);
double Create_Lagrange_Polinom(double x, std::vector<std::vector<float>> TabulPolinom);
std::vector<std::vector<float>> Lagrange_interpolate(std::vector<std::vector<float>> TabulPolinom, std::vector<float> Points);
double Create_Base_Polinom_Lagrange(double x, int indexCur, std::vector<std::vector<float>> TabulPolinom);
std::vector<std::vector<float>> Tabulate(std::vector<float> Row_of_points_X);
std::vector<float> Get_Nodes_Cheb(float beginPoint, float endPoint, int Amount_of_points);
std::vector<float> Get_Regular_Row(float beginPoint, float endPoint, int Amount_of_points);
float static Get_Step(float beginPoint, float endPoint, int Amount_of_points);
std::vector<float> Get_Info() ;
void Write_P(std::vector<std::vector<float>> dt);
void Write_P1(std::vector<float> dt);
void Write_P2(std::vector<std::vector<float>> dt);