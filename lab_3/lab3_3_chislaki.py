from lab_1 import lab1_chislaki as gauss
import matplotlib.pyplot as plt
import numpy as np


def calculate_system(matrix, size):
    with open("test.txt", "w", encoding="utf-8") as file:
        file.write(str(size) + "\n")
        for i in range(size):
            for j in range(size + 1):
                file.write(str(matrix[i][j]) + " ")
            file.write("\n")

    size, matrix, b = gauss.read_matrix("test.txt")
    P, L, U, b, permutation_counter = gauss.lu_method_forward(matrix, size, b)
    answers = gauss.lu_method_backward_vector(size, U, b)
    return answers


def mnk_method(x, y, n):
    sum_x = sum(x)
    sum_y = sum(y)
    sum_x_2 = sum(x_i ** 2 for x_i in x)
    sum_x_3 = sum(x_i ** 3 for x_i in x)
    sum_x_4 = sum(x_i ** 4 for x_i in x)
    sum_x_5 = sum(x_i ** 5 for x_i in x)
    sum_x_6 = sum(x_i ** 6 for x_i in x)
    sum_xy = sum(x_i * y_i for x_i, y_i in zip(x, y))
    sum_y_x_2 = sum(y_i * x_i ** 2 for x_i, y_i in zip(x, y))
    sum_y_x_3 = sum(y_i * x_i ** 3 for x_i, y_i in zip(x, y))

    # находим приближающий многочлен 1-ой степени
    matrix_1 =[[n, sum_x, sum_y],
               [sum_x, sum_x_2, sum_xy]]

    a_1 = calculate_system(matrix_1, 2)

    f_1 = [0 for i in range(n)]
    for i in range(n):
        f_1[i] = a_1[0] + x[i] * a_1[1]

    error_1 = 0
    for i in range(n):
        error_1 += (f_1[i] - y[i]) ** 2

    # находим приближающий многочлен 2-ой степени
    matrix_2 = [[n, sum_x, sum_x_2, sum_y],
                [sum_x, sum_x_2, sum_x_3, sum_xy],
                [sum_x_2, sum_x_3, sum_x_4, sum_y_x_2]]

    a_2 = calculate_system(matrix_2, 3)

    f_2 = [0 for i in range(n)]
    for i in range(n):
        f_2[i] = a_2[0] + x[i] * a_2[1] + a_2[2] * x[i] ** 2

    error_2 = 0
    for i in range(n):
        error_2 += (f_2[i] - y[i]) ** 2

    # находим приближающий многочлен 3-ей степени
    matrix_3 = [[n, sum_x, sum_x_2, sum_x_3, sum_y],
                [sum_x, sum_x_2, sum_x_3, sum_x_4, sum_xy],
                [sum_x_2, sum_x_3, sum_x_4, sum_x_5, sum_y_x_2],
                [sum_x_3, sum_x_4, sum_x_5, sum_x_6, sum_y_x_3]]

    a_3 = calculate_system(matrix_3, 4)

    f_3 = [0 for i in range(n)]
    for i in range(n):
        f_3[i] = a_3[0] + x[i] * a_3[1] + a_3[2] * x[i] ** 2 + a_3[3] * x[i] ** 3

    error_3 = 0
    for i in range(n):
        error_3 += (f_3[i] - y[i]) ** 2

    return a_1, a_2, a_3, f_1, f_2, f_3, error_1, error_2, error_3


def print_answer(x, a_1, a_2, a_3,
                 f_1, f_2, f_3, error_1, error_2, error_3, n):
    print(f"\nкоэффициенты многочлена 1-ой степени:\na_0 = {a_1[0]}, a_1 = {a_1[1]}\n")
    print(f"коэффициенты многочлена 2-ой степени:\na_0 = {a_2[0]}, a_1 = {a_2[1]}, a_2 = {a_2[2]}\n")
    print(f"коэффициенты многочлена 3-ей степени:\na_0 = {a_3[0]}, a_1 = {a_3[1]}, a_2 = {a_3[2]}, a_3 = {a_3[3]}\n")

    print(f"полученный многочлен P1(x):\n{a_1[0]} + {a_1[1]} * x\n")
    print(f"полученный многочлен P2(x):\n{a_2[0]} + {a_2[1]} * x + {a_2[2]} * x ** 2\n")
    print(f"полученный многочлен P3(x):\n{a_3[0]} + {a_3[1]} * x + {a_3[2]} * x ** 2 + {a_3[3]} * x ** 3\n")

    mse_1 = error_1 / n
    mse_2 = error_2 / n
    mse_3 = error_3 / n
    print(f"MSE_1 = {mse_1}, MSE_2 = {mse_2}, MSE_3 = {mse_3}\n")

    print("значения P(x_j):")
    for i in range(6):
        print("x_j =", x[i], ":", f_1[i], f_2[i], f_3[i])

    errors = [float(mse_1), float(mse_2), float(mse_3)]
    min_err = min(errors)
    num_err = 0
    for i in range(3):
        if errors[i] == min_err:
            num_err = i + 1
    print(f"\nсреднеквадратические ошибки:\n{errors}\n\nнаименьшая ошибка: {min(errors)}, многочлен степени {num_err}")


def graphic(x, y, a_1, a_2, a_3):
    x_smooth = np.linspace(min(x), max(x), 100)

    y_1_smooth = a_1[0] + a_1[1] * x_smooth
    y_2_smooth = a_2[0] + a_2[1] * x_smooth + a_2[2] * x_smooth ** 2
    y_3_smooth = a_3[0] + a_3[1] * x_smooth + a_3[2] * x_smooth ** 2 + a_3[3] * x_smooth ** 3

    plt.figure(figsize=(8, 6))

    plt.scatter(x, y, color='black', label='исходные точки', zorder=5)

    plt.plot(x_smooth, y_1_smooth, label='P1(x)', color='blue', linestyle='--')
    plt.plot(x_smooth, y_2_smooth, label='P2(x)', color='green', linestyle='-.')
    plt.plot(x_smooth, y_3_smooth, label='P3(x)', color='red', linestyle='-')

    plt.title('метод наименьших квадратов')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()


x = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
y = [-2.9502, -1.8647, -0.63212, 1.0, 3.7183, 9.3891]
n = 6
a_1, a_2, a_3, f_1, f_2, f_3, error_1, error_2, error_3 = mnk_method(x, y, n)
print_answer(x, a_1, a_2, a_3, f_1, f_2, f_3, error_1, error_2, error_3, n)

graphic(x, y, a_1, a_2, a_3)