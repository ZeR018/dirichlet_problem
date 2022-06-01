#pragma once
#include "MMN.h"
#include <iostream>
#include <omp.h>
#include <iomanip>

MMN::MMN()
{
	pi = 3.14159265359;
}


MMN::MMN(int N, int M, TVector<double> XBorder, TVector<double> YBorder, string _taskName, string _interpolation)
{
	taskName = _taskName;	// "T" или "M"
	pi = 2 * asin(1);
	n = N; m = M;

	xBorder = XBorder; yBorder = YBorder;

	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0 * n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0 * m);
	hE = -1 / pow(h, 2);
	kE = -1 / pow(k, 2);
	A = -2 * (hE + kE);

	V = TMatrix<double>(m + 1);
	U = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
	{
		V[j] = TVector<double>(n + 1);
		U[j] = TVector<double>(n + 1);
		for (int i = 0; i <= n; i++)
		{
			V[j][i] = 0.0;
			U[j][i] = 0.0;
		}
	}
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
	
	if (_interpolation == "X")
		XInterpolation();
	else if (_interpolation == "Y")
		YInterpolation();
}


double MMN::F_Function(double x, double y)
{
	double res;

	if (taskName == "M")
		res = fabs(pow(sin(pi * x * y), 3));
	else
		res = -4 * exp(1 - pow(x, 2) - pow(y, 2)) * (1 - pow(x, 2) - pow(y, 2));

	return res;
}


void MMN::FunctionInicialisation()
{
	double x, y = yBorder[0] + k;
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			F[j][i] = F_Function(x, y);
			x += h;
		}
		y += k;
	}
}


double MMN::ExactSolution(double x, double y)
{
	double res;

	res = exp(1 - pow(x, 2) - pow(y, 2));

	return res;
}


double MMN::XInicialConditions(double x, int Num)
{
	double res = 0.0;
	if (taskName == "M")
	{
		switch (Num)
		{
		case 1:
			res = fabs(sin(pi * x));
			break;
		case 2:
			res = fabs(sin(pi * x));
			break;
		}
	}
	else
	{
		switch (Num)
		{
		case 1:
			res = exp(-pow(x, 2));
			break;
		case 2:
			res = exp(-pow(x, 2));
			break;
		}
	}
	return res;
}


double MMN::YInicialConditions(double y, int Num)
{
	double res = 0.0;
	if (taskName == "M")
	{
		switch (Num)
		{
		case 1:
			res = -pow(y, 2) + 1.0;
			break;
		case 2:
			res = -pow(y, 2) + 1.0;
			break;
		}
	}
	else
	{
		switch (Num)
		{
		case 1:
			res = exp(-pow(y, 2));
			break;
		case 2:
			res = exp(-pow(y, 2));
			break;
		}
	}
	return res;
}


void MMN::Inicialisation()
{
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		V[j][0] = YInicialConditions(y, 1);
		V[j][n] = YInicialConditions(y, 2);
		U[j][0] = YInicialConditions(y, 1);
		U[j][n] = YInicialConditions(y, 2);
		y += k;
	}
	x = xBorder[0] + h;
	for (int i = 1; i < n; i++)
	{
		V[0][i] = XInicialConditions(x, 1);
		V[m][i] = XInicialConditions(x, 2);
		U[0][i] = XInicialConditions(x, 1);
		U[m][i] = XInicialConditions(x, 2);
		x += h;
	}
}


TVector<double> MMN::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name) //Точность решения
{
	return 0;
}


TVector<double> MMN::MethodError(double eps, int MaxIterations) //Погрешность решения
{
	return 0;
}


void MMN::VectorNevyazki()
{
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]) + F[j][i];
		}
}

double MMN::NevyazkaInf()
{
	double res = 0, temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = fabs(R[j][i]);
			if (temp > res)
				res = temp;
		}

	return res;
}

double MMN::NevyazkaEvkl()
{
	double res = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res += pow(R[j][i], 2);
		}

	return sqrt(res);
}

void MMN::XInterpolation()
{
	double x, a, b;
	a = xBorder[0]; b = xBorder[xBorder.Size() - 1];
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			V[j][i] = ((x - a) * V[j][n] - (x - b) * V[j][0]) / (b - a);
			x += h;
		}
	}
}

void MMN::YInterpolation()
{
	double y, c, d;
	c = yBorder[0]; d = yBorder[yBorder.Size() - 1];
	for (int i = 1; i < n; i++)
	{
		y = yBorder[0] + k;
		for (int j = 1; j < m; j++)
		{
			V[j][i] = ((y - c) * V[m][i] - (y - d) * V[0][i]) / (d - c);
			y += k;
		}
	}
}


void MMN::SaveGrid(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int i = 0; i < xBorder.Size(); i++)
		file << xBorder[i] << endl;
	file << h << endl;

	for (int i = 0; i < yBorder.Size(); i++)
		file << yBorder[i] << endl;
	file << k << endl;

	if (taskName == "T")
		file << 1 << endl;
	else
		file << 0 << endl;

	file.close();

}

void MMN::SaveData(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int j = 0; j <= m; j++)
	{

		for (int i = 0; i <= n; i++)
		{
			file << V[j][i] << endl;//"\t";
		}
		//std::cout << endl;
	}
	file.close();

}

void MMN::SaveExData(string s)
{
	ofstream file(s);
	double res;

	for (int j = 0; j <= m; j++)
	{
		for (int i = 0; i <= n; i++)
		{
			file << U[j][i] << endl;//"\t";
		}
		//std::cout << endl;
	}
	file.close();
}
//------------------------------------------------------------

void MMN::SetParams()
{
	double Arr = 0, ArAr = 0;
	double temp;
	for (int j = 1; j < m; j++)
	{
		for (int i = 1; i < n; i++)
		{
			temp = A * R[j][i] + hE * (R[j][i + 1] + R[j][i - 1]) + kE * (R[j + 1][i] + R[j - 1][i]);
			Arr += temp * R[j][i];
			ArAr += pow(temp, 2);
		}
	}
	tau = Arr / ArAr;
}

double MMN::Runner()
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];
		}
			
	SetParams();
	double accurancy = 0;

	/*int flow = 12;
	omp_set_num_threads(flow);
	#pragma omp parallel for*/
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] - tau * R[j][i];
			
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	
	return accurancy;
}

TVector<double> MMN::SolvingMainTask(TVector<double> eps, TVector<int> MaxIterations)
{
	SaveData("StartInterpolationMain.txt");
	VectorNevyazki();
	double r0_inf = NevyazkaInf();
	cout << "Start solving" << endl << endl;



	// Для половинного шага
	MMN Solution(2 * n, 2 * m, xBorder, yBorder, "M");
	Solution.SaveData("StartInterpolationMain2.txt");
	Solution.VectorNevyazki();
	double r0_inf_sup = Solution.NevyazkaInf();

	TVector<double> accurancy(2); //Точность метода
	accurancy[0] = 1.0 + eps[0];
	accurancy[1] = 1.0 + eps[1];

	TVector<int> IterationsCount(2);

	// Для обычного шага
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		//if (IterationsCount[0] % 100 == 0)
			//cout << IterationsCount[0] << "\t";
		accurancy[0] = Runner();
		IterationsCount[0]++;
	}
	cout << endl;
	cout << "exit1" << endl;

	// Для половинного шага
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		if (IterationsCount[1] % 1000 == 0)
			cout << IterationsCount[1] << "\t";
		accurancy[1] = Solution.Runner();
		IterationsCount[1]++;
	}
	cout << endl;
	cout << "exit2" << endl;

	SaveData("MainSolutA.txt"); //Заполнили файл численным решением
	Solution.SaveData("SupSolutA.txt");  //Заполнили файл с численным решением на вспомогательной сетке
	SaveGrid("DifferenceA.txt"); //Заполняем файл начальной информацией
	SaveExData("ExData.txt"); //Заполняем файл с истинным решением

	ofstream Difference("Diff.txt"); //

	double error = 0; //Точность решения
	double sup; //Вспомогательная переменная
	int ix = 0, jy = 0;
	for (int j = 0; j <= m; j++)
	{
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - Solution.V[2 * j][2 * i]);
			Difference << sup << endl;
			if (sup > error)
			{
				//cout << i << "\t" << j << "\t" << V[j][i] << "\t" << Solution.V[j][i] << endl;
				error = sup;
				ix = i;
				jy = j;
			}
		}
	}

	//for (int i = n; i >= 0; i--)
	//{
	//	for (int j = m; j >= 0; j--)
	//		cout << Solution.V[2 * j][2 * i] << '\t';
	//	cout << endl;
	//}

	double temp;
	VectorNevyazki();
	Solution.VectorNevyazki();

	sup = NevyazkaInf();
	temp = Solution.NevyazkaInf();
	TVector<double> result(11);
	result[0] = accurancy[0];			//Точность метода на сетке (n+1,m+1)
	result[1] = IterationsCount[0];		//Количество итераций на сетке (n+1,m+1)
	result[2] = accurancy[1];			//Точность метода на сетке (2n+1,2m+1)
	result[3] = IterationsCount[1];		//Количество итераций на сетке (2n+1,2m+1)
	result[4] = error;					//Точность решения
	result[5] = sup;					//Евклидова норма невязки на основной сетке
	result[6] = temp;					//Евклидова норма невязки на вспомогательной сетке
	result[7] = xBorder[0] + ix * h;	//Значение x в самой плохой точке
	result[8] = yBorder[0] + jy * k;	//Значение y в самой плохой точке
	result[9] = r0_inf;					//Невязка на начальном приближении на основной сетке
	result[10] = r0_inf_sup;			//Невязка на начальном приближении на вспомогательной сетке


	cout << "################################################# Начальные параметры #################################################" << endl;
	cout << endl;
	cout << "Программу написал Варварин Евгений, 381903_2" << endl;
	cout << "Для решения задачи на интервале по x [" << xBorder[0] << ", " << xBorder[1] << "], по y [" << yBorder[0] << ", " << yBorder[1] << "]" << endl;
	cout << "Использована равномерная сетка с числом разбиений по x n = " << n << " и числом разбиений по y m = " << m << "." << endl;
	cout << "Также была решена вспомогательная задача с числом разбиений (2n, 2m)." << endl;
	cout << "Метод: Метод Минимальных Невязок" << endl; // , параметр тау : " << tau << " (Получен на последней итерации)." << endl;
	cout << "Критерии остановки по точности = " << eps[0] << " и по числу итераций = " << MaxIterations[0] << "." << endl;
	cout << endl;
	cout << "####################################################### Справка #######################################################" << endl << endl;
	cout << "Количество итераций на основной сетке S: " << result[1] << ",  достигнутая точность: " << result[0] << "." << endl;
	cout << "Количество итераций на вспомогательной сетке S2: " << result[3] << ",  достигнутая точность: " << result[2] << "." << endl << endl;

	cout << "Для невязки СЛАУ использована норма max" << endl;

	cout << "Невязка на начальном приближении на основной сетке: " << result[9] << "." << endl;
	cout << "Невязка на начальном приближении на вспомигательной сетке: " << result[10] << "." << endl;
	cout << "На основной сетке схема решена с невязкой: " << result[5] << "." << endl;
	cout << "На вспомогательной сетке схема решена с невязкой: " << result[6] << "." << endl << endl;

	cout << "Точность решения основной задачи: " << result[4] << "." << endl;
	cout << "Максимальное отклонение между решениями на обычной и двойной сетке - x = " << result[7] << " y = " << result[8] << endl;
	cout << "Начальное приближение нулевое" << endl;
	cout << "#######################################################################################################################";
	cout << endl << endl;
	return result;
}

TVector<double> MMN::SolvingTestTask(TVector<double> eps, TVector<int> MaxIterations)
{
	SaveData("StartInterpolationTest.txt"); //Сохраняем начальное приближение
	VectorNevyazki();
	double r0_inf = NevyazkaInf();

	TVector<double> accurancy(2); //Точность метода
	accurancy[0] = 1 + eps[0];

	TVector<int> IterationsCount(2);

	// Для обычного шага
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		if (IterationsCount[0] % 1000 == 0)
			cout << IterationsCount[0] << "\t";
		accurancy[0] = Runner();
		IterationsCount[0]++;
	}
	cout << endl;
	cout << "exit" << endl;

	// Точное решение
	double y0 = yBorder[0];
	for (int j = 1; j < m; j++)
	{
		double x0 = xBorder[0];
		y0 += k;
		for (int i = 1; i < n; i++)
		{
			x0 += h;
			U[j][i] = ExactSolution(x0, y0);
		}
	}

	SaveData("MainSolutA.txt"); //Заполнили файл численным решением
	//Solution.SaveData("SupSolutA.txt");  //Заполнили файл с численным решением на вспомогательной сетке
	SaveGrid("DifferenceA.txt"); //Заполняем файл начальной информацией
	SaveExData("ExData.txt"); //Заполняем файл с истинным решением

	ofstream Difference("Diff.txt"); //

	double error = 0; //Точность решения
	double sup; //Вспомогательная переменная
	int ix = 0, jy = 0;
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - U[j][i]);
			Difference << sup << endl;
			if (sup > error)
			{
				error = sup;
				ix = i;
				jy = j;
			}
		}

	double temp;
	VectorNevyazki();
	//Solution.VectorNevyazki();

	sup = NevyazkaInf();
	//temp = Solution.NevyazkaEvkl();
	temp = 0;
	TVector<double> result(9);
	result[0] = accurancy[0];			//Точность метода на сетке (n+1,m+1)
	result[1] = IterationsCount[0];		//Количество итераций на сетке (n+1,m+1)
	result[4] = error;					//Точность решения
	result[5] = sup;					//Евклидова норма невязки на последнем шаге
	result[6] = r0_inf;				//Евклидова норма невязки на шаге 0
	result[7] = xBorder[0] + ix * h;	//Значение x в самой плохой точке
	result[8] = yBorder[0] + jy * k;	//Значение y в самой плохой точке

	cout << "################################################# Начальные параметры #################################################" << endl;
	cout << endl;
	cout << "Программу написал Варварин Евгений, 381903_2" << endl;
	cout << "Для решения задачи на интервале по x [" << xBorder[0] << ", " << xBorder[1] << "], по y [" << yBorder[0] << ", " << yBorder[1] << "]" << endl;
	cout << "Использована равномерная сетка с числом разбиений по x n = " << n << " и числом разбиений по y m = " << m << "." << endl;
	cout << "Метод: Метод Минимальных Невязок" << endl; // , параметр тау : " << tau << " (Получен на последней итерации)." << endl;
	cout << "Критерии остановки по точности = " << eps[0] << " и по числу итераций = " << MaxIterations[0] << "." << endl;
	cout << endl;
	cout << "####################################################### Справка #######################################################" << endl << endl;

	cout << "Количество итераций S: " << result[1] << ",  достигнутая точность: " << result[0] << "." << endl;
	cout << "Схема решена с невязкой: " << result[5] << " (для невязки СЛАУ использована норма max)." << endl;
	cout << "Невязка на нулевом шаге: " << result[6] << " (для невязки СЛАУ использована норма max)." << endl;
	cout << "Начальное приближение нулевое" << endl;
	cout << endl;
	cout << "Погрешность решения тестовой задачи: " << result[4] << "." << endl;
	cout << "Максимальное отклонение точного и численного решений наблюдается в узле x = " << result[7] << " y = " << result[8] << endl;
	cout << "Начальное приближение нулевое" << endl;
	cout << endl;
	cout << "#######################################################################################################################";
	cout << endl << endl;
	return result;
}

