#include <iostream>
#include "TMatrix.h"
#include "MMN.h"
#include <Windows.h>

int main()
{
    setlocale(LC_ALL, "Russian");

    // Параметры 
    int N = 40;
    int M = 40;
    int Max_iter = 10000000;          // Максимальное число шагов
    double er = 1e-11;            // Остановка поточности
    TVector <double> xBorder(2);
    TVector <double> yBorder(2);
    TVector <double> accurancy(2);
    TVector <int> max_iter(2);
    accurancy[0] = er;
    accurancy[1] = er/10.0;
    max_iter[0] = Max_iter;
    max_iter[1] = Max_iter;
    xBorder[0] = -1;                // a
    xBorder[1] = 1;                 // b 
    yBorder[0] = -1;                // c
    yBorder[1] = 1;                 // d
    
    int k;
    cout << "Выберите задачу" << endl;
    cout << "Основная задача - 1" << endl;
    cout << "Тестовая задача - 2," << endl;
    cin >> k;

    // Main task
    if (k == 1)
    {
        MMN mat(N, M, xBorder, yBorder, "M");
        mat.SolvingMainTask(accurancy, max_iter);
        system("python graph.py");
    }
    // Test task
    else if (k == 2)
    {
        MMN mat(N, M, xBorder, yBorder, "T");
        mat.SolvingTestTask(accurancy, max_iter);
        system("python graph.py");
    }
    else
        cout << "Некорректный параметр. Перезапустите программу";

}