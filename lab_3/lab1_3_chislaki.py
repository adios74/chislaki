import math
import matplotlib.pyplot as plt
import numpy as np


def calculate_func(x):
    return math.exp(x) + x


def newton_polinom(x, y, x_p):
    a = [0.0 for _ in range(4)]

    a[0] = y[0]
    a[1] = (y[1] - a[0]) / (x[1] - x[0])
    a[2] = (y[2] - a[0] - a[1] * (x[2] - x[0])) / ((x[2] - x[0]) * (x[2] - x[1]))
    a[3] = (y[3] - a[0] - a[1] * (x[3] - x[0]) - a[2] * (x[3] - x[0]) * (x[3] - x[1])) / ((x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2]))

    table = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        table[i][0] = y[i]

    for j in range(1, 4):
        for i in range(4 - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i])


    polinom = (f"{a[0]} + {a[1]} * (x - {x[0]}) \n"
               f"+ {a[2]} * (x - {x[0]}) * (x - {x[1]}) + {a[3]} * (x - {x[0]}) * (x - {x[1]}) * (x - {x[2]})")

    func_x_precise = a[0] + a[1] * (x_p - x[0]) + a[2] * (x_p - x[0]) * (x_p - x[1]) + a[3] * (x_p - x[0]) * (x_p - x[1]) * (x_p - x[2])

    test = calculate_func(x_p)
    delta = abs(test - func_x_precise)
    return func_x_precise, test, delta, a, polinom, table


def lagrange_polinom(x, y, x_p):
    a = [0.0 for _ in range(4)]

    znamen = []

    znamen.append((x[0] - x[1]) * (x[0] - x[2]) * (x[0] - x[3]))
    znamen.append((x[1] - x[0]) * (x[1] - x[2]) * (x[1] - x[3]))
    znamen.append((x[2] - x[0]) * (x[2] - x[1]) * (x[2] - x[3]))
    znamen.append((x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2]))

    for i in range(4):
        a[i] = y[i] / znamen[i]

    polinom = (f"{a[0]} * (x - {x[1]}) * (x - {x[2]}) * (x - {x[3]}) +\n"
                f"+ {a[1]} * (x - {x[0]}) * (x - {x[2]}) * (x - {x[3]}) +\n"
                f"+ {a[2]} * (x - {x[0]}) * (x - {x[1]}) * (x - {x[3]}) +\n"
                f"+ {a[3]} * (x - {x[0]}) * (x - {x[1]}) * (x - {x[2]})\n")

    func_x_precise = a[0] * (x_p - x[1]) * (x_p - x[2]) * (x_p - x[3])\
    + a[1] * (x_p - x[0]) * (x_p - x[2]) * (x_p - x[3])\
    + a[2] * (x_p - x[0]) * (x_p - x[1]) * (x_p - x[3])\
    + a[3] * (x_p - x[0]) * (x_p - x[1]) * (x_p - x[2])

    test = calculate_func(x_p)
    delta = abs(test - func_x_precise)

    return func_x_precise, test, delta, a, polinom


def print_answer(x, x_p):
    y = [calculate_func(x) for x in x]
    print(f"исходные узлы интерполяции: {x}\n")
    print(f"значения функции y: {y}\n")
    func_x_precise_l, test_l, delta_l, a_l, polinom_l = lagrange_polinom(x, y, x_p)
    func_x_precise_n, test_n, delta_n, a_n, polinom_n, table = newton_polinom(x, y, x_p)
    print(f"интерполяционный многочлен Лагранжа:\n{polinom_l}")

    print("таблица разделенных разностей:")
    for i in range(len(x)):
        row = [f"{x[i]:.2f}"] + [f"{table[i][j]:.6f}" for j in range(len(x) - i)]
        print(*row, sep='\t')

    print(f"\nинтерполяционный многочлен Ньютона:\n{polinom_n}\n")

    print(f"значение многочлена Лагранжа в заданной точке x* = -0.5:\n{func_x_precise_l}\n")
    print(f"значение многочлена Ньютона в заданной точке x* = -0.5:\n{func_x_precise_n}\n")

    print(f"точное значение функции f(x*):\n{test_l}\n")
    print(f"абсолютная погрешность интерполяции:\n{delta_l}\n{delta_n}\n")

    print(f"разность между результатами методов Лагранжа и Ньютона:\n{abs(func_x_precise_l - func_x_precise_n)}\n")


def graphic(x, x_p):
    y = [calculate_func(xi) for xi in x]

    _, _, _, a_l, _ = lagrange_polinom(x, y, x_p)
    _, _, _, a_n, _, _ = newton_polinom(x, y, x_p)

    x_smooth = np.linspace(min(x) - 0.5, max(x) + 0.5, 200)

    y_func = [calculate_func(xi) for xi in x_smooth]

    y_newton = []
    for xi in x_smooth:
        val = (a_n[0] + a_n[1] * (xi - x[0]) +
               a_n[2] * (xi - x[0]) * (xi - x[1]) +
               a_n[3] * (xi - x[0]) * (xi - x[1]) * (xi - x[2]))
        y_newton.append(val)

    y_lagrange = []
    for xi in x_smooth:
        val = (a_l[0] * (xi - x[1]) * (xi - x[2]) * (xi - x[3]) +
               a_l[1] * (xi - x[0]) * (xi - x[2]) * (xi - x[3]) +
               a_l[2] * (xi - x[0]) * (xi - x[1]) * (xi - x[3]) +
               a_l[3] * (xi - x[0]) * (xi - x[1]) * (xi - x[2]))
        y_lagrange.append(val)

    plt.figure(figsize=(10, 6))

    plt.plot(x_smooth, y_func, color='black', label='исходная функция f(x)', linewidth=2)
    plt.plot(x_smooth, y_lagrange, color='red', linestyle='--', label='многочлен Лагранжа', linewidth=2)
    plt.plot(x_smooth, y_newton, color='blue', linestyle=':', label='многочлен Ньютона', linewidth=3)

    plt.scatter(x, y, color='black', zorder=5, label='узлы интерполяции')

    y_xp = calculate_func(x_p)
    plt.scatter(x_p, y_xp, color='green', label=f'точка x* = {x_p}')

    plt.title('интерполяция функции (Лагранж и Ньютон)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()


x_a = [-2, -1, 0, 1]
x_b = [-2, -1, 0.2, 1]
print("==================================================\n")
print("результат выполнения пункта а):\n")
print_answer(x_a, -0.5)
graphic(x_a, -0.5)
print("==================================================\n")
print("результат выполнения пункта б):\n")
print_answer(x_b, -0.5)
graphic(x_b, -0.5)

