import matplotlib.pyplot as plt
import pylab
#from prettytable import PrettyTable
import numpy as np
import pandas as pd
import openpyxl


tr_y = []
tr_dif = []
tr_ex = []
tr_sup = []

# Шарим файлы
with open("MainSolutA.txt") as f:
    for line in f:
        tr_y.append([float(x) for x in line.split()])

with open("DifferenceA.txt") as d:
    for line in d:
        tr_dif.append([float(x) for x in line.split()])

with open("ExData.txt") as e:
    for line in e:
        tr_ex.append([float(x) for x in line.split()])

with open("SupSolutA.txt") as e:
    for line in e:
        tr_sup.append([float(x) for x in line.split()])

x = np.array([])
xs = np.array([])
xs2 = np.array([])
y = np.array([])
th = []
z = np.array([])
ex = np.array([])
ex_data = np.array([])

# Записываем n,m DifferenceA
n_ = tr_y[0][0]
m_ = tr_y[1][0]
n = int(n_)
m = int(m_)

# численное решение z на сетке
for i in range(2, len(tr_y)):
    z = np.append(z, tr_y[i])

for i in range(0, len(tr_ex)):
    ex = np.append(ex, tr_ex[i])

# Записываем параметры из файла DifferenceA
h = tr_dif[4][0]
k = tr_dif[7][0]
x_border1 = tr_dif[2][0]
x_border2 = tr_dif[3][0]
y_border1 = tr_dif[5][0]
y_border2 = tr_dif[6][0]
taskName = tr_dif[8][0]     # 1 - test task, 0 - main task

# Составляем массивы для графиков на основной сетке
x_ = x_border1
while x_ <= x_border2:
    xs = np.append(xs, round(x_, 5))
    x_ += h
    x_ = round(x_, 5)

for i in range(n+1):
    j = 0
    while j < (m+1):
        y = np.append(y, round(y_border1 + j * h, 5))
        j+=1

for j in range(len(xs)):
    for i in range(n+1):
        x = np.append(x, xs[j])

# Test task
if taskName == 1:
    # точное решение u* (x, y)(поверхность);
    # начальное приближение v(0) (xi, yj)(поверхность);
    # численное решение v(N) (x, y)(поверхность);
    # разность точного и численного решения(поверхность).

    # Разность решений
    tr_difference_eq = []
    with open("Diff.txt") as f:
        for line in f:
            tr_difference_eq.append([float(x) for x in line.split()])
    difference = np.array([])
    for i in range(0, len(tr_difference_eq)):
        difference = np.append(difference, tr_difference_eq[i])
    # Для таблицы разности решений
    diff_data = np.array([])
    for i in range(len(z)):
        diff_data = np.append(diff_data, difference[i])
    # Начальное приближение
    tr_s_interp = []
    with open("StartInterpolationTest.txt") as f:
        for line in f:
            tr_s_interp.append([float(x) for x in line.split()])
    s_interp = []
    for i in range(2, len(tr_s_interp)):
        s_interp.append(tr_s_interp[i][0])

    # Для таблицы численного решения
    td_data = np.array([])
    for i in range(len(z)):
        td_data = np.append(td_data, z[i])

    for i in range(len(ex)):
        ex_data = np.append(ex_data, ex[i])

    # Увеличить сетку на 2 по x и по y
    # z_int_test = np.array([]) # start interpolation array test
    # print(len(s_interp))
    # iterator = 0
    # for j in range(0, m+1):
    #     for i in range(0, n+1):
    #         if (i == 1 and j != 0 and j != m+2) or (i == n + 1 and j != 0 and j != m+2)\
    #                 or (j == 1 and i != 0 and i != n+2) or (j == m + 1 and i != 0 and i != n+2):
    #             z_int_test = np.append(z_int_test, 0)
    #             iterator+=1
    #         z_int_test = np.append(z_int_test, s_interp[i + j * n + 1])
    # print("it: ", iterator)
    #
    # x_int_test = np.array([])
    # for j in range(len(xs)):
    #     x_int_test = np.append(x_int_test, x_border1)
    #     for i in range(n + 1):
    #         x_int_test = np.append(x_int_test, xs[j])
    #     x_int_test = np.append(x_int_test, x_border2)
    #
    # y_int_test = np.array([])
    # for i in range(n + 1):
    #     j = 0
    #     while j < (m + 1):
    #         y_int_test = np.append(y_int_test, round(y_border1 + j * h, 5))
    #         if (j == 0):
    #             y_int_test = np.append(y_int_test, round(y_border1, 5))
    #         elif (j == m):
    #             y_int_test = np.append(y_int_test, round(y_border2, 5))
    #         j+=1
    # print(len(x_int_test), len(y_int_test), len(z_int_test))

    ####################### Таблицы
    # Точное решение
    ex_data = ex_data.reshape(int(m) + 1, int(n) + 1)
    df_ex = pd.DataFrame(ex_data)
    df_ex.to_excel("Table_ex.xlsx")
    # Численное решение
    td_data = td_data.reshape(int(m) + 1, int(n) + 1)
    df = pd.DataFrame(td_data)
    df.to_excel("Table_chisl.xlsx")
    # Разность точного решения и численного
    diff_data = diff_data.reshape(int(m) + 1, int(n) + 1)
    diff_data = pd.DataFrame(diff_data)
    diff_data.to_excel("Table_diff.xlsx")

    ####################### Графики
    # Истинное
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(x, y, ex)
    ax.set_title('Точное решение')
    # Численное
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.plot_trisurf(x, y, z)
    ax2.set_title('Численное решение')
    # Разность истинного и численного
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.plot_trisurf(x, y, difference)
    ax3.set_title('Разность точного решения и численного')
    # Начальное приближение
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111, projection='3d')
    ax4.plot_trisurf(x, y, s_interp)
    ax4.set_title('Начальное приближение')

    plt.show()
elif taskName == 0:
    # начальное приближение v(0)(xi, yj)(поверхность);
    # начальное приближение v2(0)(xi, yj)(поверхность);
    # численное решение v(N) (x, y)(поверхность);
    # численное решение v2(N2)(x, y), полученное на сетке с половинным шагом (поверхность);
    # разность численных решений v(N) (x, y) и v2(N2)(x, y)(поверхность).

    # Записываем решение z вспомогательной задачи
    z_sup = []
    n2 = (int)(tr_sup[0][0])
    m2 = (int)(tr_sup[1][0])
    h2 = h / 2.0
    k2 = k / 2.0
    for i in range(2, len(tr_sup)):
        z_sup = np.append(z_sup, tr_sup[i])

    xs_sup = []
    x_sup = []
    y_sup = []

    # Составляем массивы для графиков на вспомогательной сетке
    x_ = x_border1
    while x_ <= x_border2:
        xs_sup = np.append(xs_sup, round(x_, 5))
        x_ += h2
        x_ = round(x_, 5)

    for i in range(n2 + 1):
        j = 0
        while j < (m2 + 1):
            y_sup = np.append(y_sup, round(y_border1 + j * k2, 5))
            j += 1

    for j in range(len(xs_sup)):
        for i in range(n2 + 1):
            x_sup = np.append(x_sup, xs_sup[j])

    # Разница решений на основной и вспомогательной сетках
    tr_difference_eq = []
    with open("Diff.txt") as f:
        for line in f:
            tr_difference_eq.append([float(x) for x in line.split()])
    difference = []
    for i in range(0, len(tr_difference_eq)):
        difference.append(tr_difference_eq[i][0])

    # Для таблицы численного решения
    td_data = np.array([])
    diff_data = np.array([])
    for i in range(len(z)):
        td_data = np.append(td_data, z[i])

    for i in range(len(difference)):
        diff_data = np.append(diff_data, difference[i])

    z_sup_data = np.array([])
    for i in range(len(z_sup)):
        z_sup_data = np.append(z_sup_data, z_sup[i])

    # Начальные приближения
    tr_start_interp_main_1 = []
    tr_start_interp_main_2 = []
    with open("StartInterpolationMain.txt") as f:
        for line in f:
            tr_start_interp_main_1.append([float(x) for x in line.split()])
    with open("StartInterpolationMain2.txt") as f:
        for line in f:
            tr_start_interp_main_2.append([float(x) for x in line.split()])

    s_interp_m_1 = []
    s_interp_m_2 = []
    for i in range(2, len(tr_start_interp_main_1)):
        s_interp_m_1.append(tr_start_interp_main_1[i][0])

    for i in range(2, len(tr_start_interp_main_2)):
        s_interp_m_2.append(tr_start_interp_main_2[i][0])

    ####################### Таблицы
    # Численное решение
    td_data = td_data.reshape(int(m) + 1, int(n) + 1)
    df = pd.DataFrame(td_data)
    df.to_excel("Table_chisl.xlsx")
    # Численное решение на вспомогательной сетке
    z_sup_data = z_sup_data.reshape(2 * int(m) + 1, 2 * int(n) + 1)
    z_sup_data = pd.DataFrame(z_sup_data)
    z_sup_data.to_excel("Table_chisl_2.xlsx")
    # Разность точного решения и численного
    diff_data = diff_data.reshape(int(m) + 1, int(n) + 1)
    diff_data = pd.DataFrame(diff_data)
    diff_data.to_excel("Table_diff.xlsx")

    ####################### Графики

    # Численное
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(x, y, z)
    ax.set_title('Численное решение основной задачи')
    # Численное на вспомогательной сетке
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.plot_trisurf(x_sup, y_sup, z_sup)
    ax2.set_title('Численное решение вспомогательной задачи')
    # Разница точного решения и численного
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.plot_trisurf(x, y, difference)
    ax3.set_title('Разница решений основной и вспомогательной задач')
    # Начальное приближение на основной сетке
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111, projection='3d')
    ax4.plot_trisurf(x, y, s_interp_m_1)
    ax4.set_title('Начальное приближение для основной задачи')
    # Начальное приближение на вспомогательной сетке
    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111, projection='3d')
    ax5.plot_trisurf(x_sup, y_sup, s_interp_m_2)
    ax5.set_title('Начальное приближение для дополнительной задачи')
    plt.show()
else:
    print('data error')
    exit(0)

exit(0)
