#pragma once

#ifndef MMN_H
#define MMN_H

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include "TMatrix.h"

using namespace std;

class MMN
{
protected:
	TMatrix<double> V, R, F, U;			//������� ������, ������� � �������� ������� � ������ ����� ��������������
	TVector<double> xBorder, yBorder;	//������� �� ��� x � y ��������������
	double h, k;						//��� �� ��� x � y ��������������
	double hE, kE, A;					//��������������� ������
	int n, m;							//����������� �����
	double pi;							// ��
	string taskName;					// "T" ��� "M" - �������� ������ ��� �������� ������

public:
	MMN();
	MMN(int N, int M, TVector<double> XBorder, TVector<double> YBorder, string _taskName = "T", string _interpolation = "0"); //����������� �������������

	double F_Function(double x, double y); //������� � ������ ����� (-f(x,y))
	void FunctionInicialisation();
	double ExactSolution(double x, double y);//������ �������
	double XInicialConditions(double x, int Num); //��������� ������� � ���� ������� ����������� ��� X
	double YInicialConditions(double y, int Num); //��������� ������� � ���� ������� ����������� ��� Y
	void Inicialisation(); //���������� ��������� ������� � ���� �������

	//virtual double Runner(); //������������ �����


	virtual TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name); //�������� �������, ��� �������� �����
	virtual TVector<double> MethodError(double eps, int MaxIterations); //����������� �������, ��� �������� �����

	void VectorNevyazki();
	double NevyazkaInf(); //������� �� ����� �������������
	double NevyazkaEvkl(); //��������� ����� �������

	void XInterpolation();
	void YInterpolation();

	void SaveGrid(string s);
	void SaveData(string s);
	void SaveExData(string s);

	//----------------------------------------------------------------------------------
private:
	double tau;
public:
	void SetParams();
	double Runner();
	TVector<double> SolvingMainTask(TVector<double> eps, TVector<int> MaxIterations);
	TVector<double> SolvingTestTask(TVector<double> eps, TVector<int> MaxIterations);
};

#endif
