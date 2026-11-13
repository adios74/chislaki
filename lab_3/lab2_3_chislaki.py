from lab_1 import lab2_chislaki
import matplotlib.pyplot as plt
import numpy as np


def splines_method(x, y, n, h, x_precise):
    # с помощью формул 3.13 и метода прогонки получаем коэффициенты с2, с3, с4
    tridiagonal_matrix = [[0] * n for _ in range(n - 1)]
    for i in range(n - 1):
        tridiagonal_matrix[i][i] = 2 * (h + h)
        tridiagonal_matrix[i][n - 1] = 3 * ((y[i + 2] - y[i + 1]) / h - (y[i + 1] - y[i]) / h)
        for j in range(n - 1):
            if j == i - 1 or j == i + 1:
                tridiagonal_matrix[i][j] = h

    with open("test.txt", "w", encoding="utf-8") as file:
        file.write(str(n - 1) + "\n")
        for i in range(n - 1):
            for j in range(n):
                file.write(str(tridiagonal_matrix[i][j]) + " ")
            file.write("\n")

    matrix, size = lab2_chislaki.read_elements("test.txt")
    answers, coefficients = lab2_chislaki.tridiagonal_method(matrix, size)
    c = [0 for _ in range(n + 1)]
    for i in range(n + 1):
        if i == 0 or i == n:
            c[i] = 0.0
        else:
            c[i] = answers[i - 1]

    # с помощью формул 3.14 получаем остальные коэффициенты
    a = [0 for _ in range(n)]
    b = [0 for _ in range(n)]
    d = [0 for _ in range(n)]

    for i in range(n):
        a[i] = y[i]

        if i == n - 1:
            b[i] = (y[i + 1] - y[i]) / h - 2/3 * h * c[i]
            d[i] = - c[i] / (3 * h)
        else:
            b[i] = (y[i + 1] - y[i]) / h - 1/3 * h * (c[i + 1] + 2 * c[i])
            d[i] = (c[i + 1] - c[i]) / (3 * h)

    # находим нужный отрезок и вычисляем значение для x* = -0.5
    segment = 0
    for i in range(len(x)):
        if (i != len(x) - 1):
            if x[i] <= x_precise <= x[i + 1]:
                segment = i
                break
        else:
            if x[i] <= x_precise:
                print("ошибка: x не лежит в отрезке")

    f_x_p = (a[segment] + b[segment] * (x_precise - x[segment]) + c[segment] * (x_precise - x[segment]) ** 2
             + d[segment] * (x_precise - x[segment]) ** 3)

    return a, b, c, d, segment, f_x_p


def print_answer(x, y, a, b, c, d, segment, f_x_p):
    n = 4
    print(f"исходные узлы:\n{x}\n")
    print(f"исходные значения функции:\n{y}\n")

    print(f"коэффициенты и интервалы, на которых действует соответствующий коэффициент:")
    print(f"i\t[x(i); x(i + 1)]\ta[i]\tb[i]\tc[i]\td[i]")
    for i in range(n):
        print(f"{i}\t[{x[i]}; {x[i + 1]}]\t{a[i]:.5f}\t{b[i]:.5f}\t{c[i]:.5f}\t{d[i]:.5f}")

    print(f"\nзначение сплайна в точке x* = -0.5:\n{f_x_p}\n")
    print(f"номер интервала, содержащего точку x*:\n{segment}\n")

    ok = True
    epsilon = 1e-6
    for i in range(n):
        left = a[i]
        if abs(left - y[i]) > epsilon:
            ok = False
            break

        right = a[i] + b[i] * h + c[i] * (h ** 2) + d[i] * (h ** 3)
        if abs(right - y[i + 1]) > epsilon:
            ok = False
            break

    if ok:
        print("проверка прохождения через исходные узлы пройдена!")
    else:
        print("ошибка! проверка прохождения через исходные узлы не пройдена!")

    d_l_2 = 2 * c[0] + 6 * d[0] * (x[0] - x[0])
    d_r_2 = 2 * c[n - 1] + 6 * d[n - 1] * (x[n] - x[n - 1])

    if abs(d_l_2) < epsilon and abs(d_r_2) < epsilon:
        print("проверка вторых производных пройдена!")
    else:
        print("ошибка! проверка вторых производных не пройдена!")


def graphic(x, y, a, b, c, d, x_precise, f_x_p):
    plt.figure(figsize=(8, 6))
    n = len(x) - 1
    for i in range(n):
        x_segment = np.linspace(x[i], x[i + 1], 20)

        y_segment = (a[i] + b[i] * (x_segment - x[i]) +
                     c[i] * (x_segment - x[i]) ** 2 +
                     d[i] * (x_segment - x[i]) ** 3)

        if i == 0:
            plt.plot(x_segment, y_segment, color='blue', label='кубический сплайн')
        else:
            plt.plot(x_segment, y_segment, color='blue')

    plt.scatter(x, y, color='black', zorder=5, label='исходные узлы')

    plt.scatter(x_precise, f_x_p, color='red', label=f'точка x* = {x_precise}')

    plt.title('интерполяция кубическим сплайном')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()


x = [-2.0, -1.0, 0.0, 1.0, 2.0]
y = [-1.8647, -0.63212, 1.0, 3.7183, 9.3891]
n, h = 4, 1
x_precise = -0.5
a, b, c, d, segment, f_x_p = splines_method(x, y, n, h, x_precise)
print_answer(x, y, a, b, c, d, segment, f_x_p)

graphic(x, y, a, b, c, d, x_precise, f_x_p)