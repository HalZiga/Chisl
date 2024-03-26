#include "Task.hpp"

int main()
{
	std::vector<float> data = Get_Info();
	
	Write_P2(Experiment_Lagrange_Cheb(data));
	std::cout << std::endl;
	Write_P2(Experiment_Lagrange_Row(data));
}

void Write_P2(std::vector<std::vector<float>> dt)
{
	std::string writer = "";
	for (size_t i = 0; i < dt.size(); i++)
	{
		writer += std::to_string(dt[i][0]) + "|";
	}
	writer += "\n";
	for (size_t i = 0; i < dt.size(); i++)
	{
		writer += std::to_string(dt[i][1]) + "|";
	}
	writer += "\n";
	std::cout << writer;
}

void Write_P(std::vector<std::vector<float>> dt)
{
	std::string writer = "";
	for (size_t i = 0; i < dt[0].size(); i++)
	{
		writer += std::to_string(dt[0][i]) + "|";
	}
	writer += "\n";
	for (size_t i = 0; i < dt[1].size(); i++)
	{
		writer += std::to_string(dt[1][i]) + "|";
	}
	writer += "\n";
	std::cout << writer;
}

void Write_P1(std::vector<float> dt)
{
	std::string writer = "";
	for (size_t i = 0; i < dt.size(); i++)
	{
		writer += std::to_string(dt[i]) + "|";
	}
	writer += "\n";
	std::cout << writer;
}

std::vector<float> Get_Info() 
{
	std::vector<float> data;
	float a;
	std::cout << "Точка начала = ";
	std::cin >> a;
	data.push_back(a);
	std::cout << "Точка конца = ";
	std::cin >> a;
	data.push_back(a);
	std::cout << "Количество точек при табуляции = ";
	std::cin >> a;
	data.push_back(a);
	std::cout << "Количество точек при интерполяции = ";
	std::cin >> a;
	data.push_back(a);
	std::cout << "Начальное количество точек в опыте = ";
	std::cin >> a;
	data.push_back(a);
	std::cout << "Конечное количество точек в опыте = ";
	std::cin >> a;
	data.push_back(a);
	return data;
}

float static Get_Step(float beginPoint, float endPoint, int Amount_of_points)
{
	return (endPoint - beginPoint) / Amount_of_points;
}

std::vector<float> Get_Regular_Row(float beginPoint, float endPoint, int Amount_of_points)
{
	std::vector<float> Nodes;
	float step = Get_Step(beginPoint, endPoint, Amount_of_points);
	for (float x = beginPoint; x < endPoint; x += step)
	{
		Nodes.push_back(x);
	}
	Nodes.push_back(endPoint);
	return Nodes;
}

std::vector<float> Cheb_Nods(float beginPoint, float endPoint, int Amount_of_points)
{
	std::vector<float> Nodes = Get_Regular_Row(beginPoint, endPoint, Amount_of_points);
	std::vector<float> Nodes_Cheb;
	float a = beginPoint, b = endPoint, n = Amount_of_points;
	if(beginPoint != 0)
	{
		for (int i = 1; i <= Nodes.size(); i++)
		{
			double node = (a + b) / 2.0 + (b - a) / 2.0 * cos((2.0 * i - 1.0) * M_PI / (2.0 * n));
			Nodes_Cheb.push_back(node);
		}
	}
	else
	{
		for (int i = 0; i < Nodes.size(); i++)
		{
			Nodes_Cheb.push_back(0.5 * (a + b) + 0.5 * (b - a) * cos(((2.0 * i + 1) / (2.0 * n + 2.0)) * M_PI));
		}
		reverse(Nodes_Cheb.begin(), Nodes_Cheb.end());
	}
	return Nodes_Cheb;
}

std::vector<std::vector<float>> Tabulate(std::vector<float> Row_of_points_X)
{
	std::vector<float> ValueNodes;
	std::vector<std::vector<float>> result;
	for (int i = 0; i < Row_of_points_X.size(); i++)
	{
		float sum, n, a, q;
		float x = Row_of_points_X[i];
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
		ValueNodes.push_back(sum);
	}
	result.push_back(Row_of_points_X);
	result.push_back(ValueNodes);
	return result;
}

double Create_Base_Polinom_Lagrange(double x, int indexCur, std::vector<std::vector<float>> TabulPolinom)
{
	double result = 1;
	int count = TabulPolinom[0].size();
	for (int j = 0; j < count; j++)
	{
		if (j != indexCur)
		{
			result *= (x - TabulPolinom[0][j]);
			result /= (TabulPolinom[0][indexCur] - TabulPolinom[0][j]);
		}
	}
	return result;
}

std::vector<std::vector<float>> Lagrange_interpolate(std::vector<std::vector<float>> TabulPolinom, std::vector<float> Points)
{
	std::vector<float> InterPolinomF;
	for(int i = 0; i < Points.size(); i++)
	{
		InterPolinomF.push_back(Create_Lagrange_Polinom(Points[i], TabulPolinom));
	}
	std::vector<std::vector<float>> interpolate_polinom = { Points, InterPolinomF };
	return interpolate_polinom;
}

double Create_Lagrange_Polinom(double x, std::vector<std::vector<float>> TabulPolinom)
{
	int count = TabulPolinom[0].size();
	double result = 0;
	for (int i = 0; i < count; i++)
	{
		result += TabulPolinom[1][i] * Create_Base_Polinom_Lagrange(x, i, TabulPolinom);
	}
	return result;
}

std::vector<std::vector<float>> Get_Pogreshnost(std::vector<std::vector<float>> TabulatePolinom, std::vector<std::vector<float>> InterpolatePolinom)
{
	std::vector<float> Pogreshnost;
	for (size_t i = 0; i < InterpolatePolinom[1].size(); i++)
	{
		Pogreshnost.push_back(fabs(TabulatePolinom[1][i] - InterpolatePolinom[1][i]));
	}
	std::vector<std::vector<float>> result = { TabulatePolinom[0], Pogreshnost };// не совсем
	return result;
}

float Get_Max_Pogreshnost(std::vector<std::vector<float>> Pogreshnost)
{
	float max = Pogreshnost[1][0], point = Pogreshnost[0][0], count = Pogreshnost[1].size();
	for (float i = 1; i < count; i++)
	{
		if (max < Pogreshnost[1][i])
		{
			max = Pogreshnost[1][i];
			point = Pogreshnost[0][i];
		}
	}
	return max;
}

std::vector<std::vector<float>> Function_Lagrange_Row(std::vector<float> data)
{
	std::vector<float> nodes = Get_Regular_Row(data[0], data[1], data[2]);//обычное разбиение
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	std::vector<float> points_interpolate = Get_Regular_Row(data[0], data[1], data[3]);
	std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
	return interpolate_polinom;
}

std::vector<std::vector<float>> Function_Lagrange_Cheb(std::vector<float> data)
{
	std::vector<float> nodes = Cheb_Nods(data[0], data[1], data[2]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	std::vector<float> points_interpolate = Cheb_Nods(data[0], data[1], data[3]);
	std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
	return interpolate_polinom;
}

std::vector<std::vector<float>> Function_Nodes_Row(std::vector<float> data)
{
	std::vector<float> nodes = Get_Regular_Row(data[0], data[1], data[2]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	return tabulate_polinom;
}

std::vector<std::vector<float>> Function_Nodes_Cheb(std::vector<float> data)
{
	std::vector<float> nodes = Cheb_Nods(data[0], data[1], data[2]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	return tabulate_polinom;
}

std::vector<std::vector<float>> Function_Nodes_Row_2(std::vector<float> data)
{
	std::vector<float> nodes = Get_Regular_Row(data[0], data[1], data[3]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	return tabulate_polinom;
}

std::vector<std::vector<float>> Function_Nodes_Cheb_2(std::vector<float> data)
{
	std::vector<float> nodes = Cheb_Nods(data[0], data[1], data[3]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	return tabulate_polinom;
}

std::vector<std::vector<float>> Pogreshnost_Lagrange_Cheb(std::vector<float> data)
{
	std::vector<float> nodes = Cheb_Nods(data[0], data[1], data[2]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	std::vector<float> points_interpolate = Cheb_Nods(data[0], data[1], data[3]);
	std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
	nodes = Cheb_Nods(data[0], data[1], data[3]);
	tabulate_polinom = Tabulate(nodes);
	std::vector<std::vector<float>> pogreshnost = Get_Pogreshnost(tabulate_polinom, interpolate_polinom);
	return pogreshnost;
}

std::vector<std::vector<float>> Pogreshnost_Lagrange_Row(std::vector<float> data)
{
	std::vector<float> nodes = Get_Regular_Row(data[0], data[1], data[2]);
	std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
	std::vector<float> points_interpolate = Get_Regular_Row(data[0], data[1], data[3]);
	std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
	nodes = Get_Regular_Row(data[0], data[1], data[3]);
	tabulate_polinom = Tabulate(nodes);
	std::vector<std::vector<float>> pogreshnost = Get_Pogreshnost(tabulate_polinom, interpolate_polinom);
	return pogreshnost;
}

std::vector<std::vector<float>> Experiment_Lagrange_Row(std::vector<float> data)
{
	std::vector<std::vector<float>> experiment_result;
	for (float i = data[4]; i < data[5]; i++)
	{
		std::vector<float> nodes = Get_Regular_Row(data[0], data[1], i);
		std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
		std::vector<float> points_interpolate = Get_Regular_Row(data[0], data[1], data[3]);
		std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
		nodes = Get_Regular_Row(data[0], data[1], data[3]);
		tabulate_polinom = Tabulate(nodes);
		experiment_result.push_back({ i + 1, Get_Max_Pogreshnost(Get_Pogreshnost(tabulate_polinom, interpolate_polinom)) });
	}
	return experiment_result;
}

std::vector<std::vector<float>> Experiment_Lagrange_Cheb(std::vector<float> data)
{
	std::vector<std::vector<float>> experiment_result;
	for (float i = data[4]; i < data[5]; i++)
	{
		std::vector<float> nodes = Cheb_Nods(data[0], data[1], i);//ïîëó÷èëè ðàçáèåíèå îòðåçêà
		std::vector<std::vector<float>> tabulate_polinom = Tabulate(nodes);
		std::vector<float> points_interpolate = Cheb_Nods(data[0], data[1], data[3]);
		std::vector<std::vector<float>> interpolate_polinom = Lagrange_interpolate(tabulate_polinom, points_interpolate);
		nodes = Cheb_Nods(data[0], data[1], data[3]);
		tabulate_polinom = Tabulate(nodes);
		experiment_result.push_back({ i + 1, Get_Max_Pogreshnost(Get_Pogreshnost(tabulate_polinom, interpolate_polinom)) });
	}
	return experiment_result;
}
