import math
import matplotlib.pyplot as plt


def runge_romberg_method(y_h, y_2h, p):
    return (y_h - y_2h) / (2 ** p - 1)


def func_y(y, z, h):
    return z


def func_z(x, y, z, h):
    return ((x + 1) * z - y) / x


def euler_method(x_start, x_end, h):
    steps = int(round((x_end - x_start) / h))

    y = [0] * (steps + 1)
    z = [0] * (steps + 1)
    x = [x_start + i * h for i in range(steps + 1)]
    precise_ans = [0] * (steps + 1)

    # заменяем уравнение 2-го порядка системой ур-ний
    y[0] = 2 + math.exp(1)
    z[0] = 1 + math.exp(1)

    # y(i) = y(i - 1) + h * f(x(i - 1), y(i - 1), z(i - 1))
    # z(i) = z(i - 1) + h * g(x(i - 1), y(i - 1), z(i - 1))
    for i in range(steps):
        y[i + 1] = y[i] + h * func_y(y[i], z[i], h)
        z[i + 1] = z[i] + h * func_z(x[i], y[i], z[i], h)
        precise_ans[i] = x[i] + 1 + math.exp(x[i])

    precise_ans[-1] = x[-1] + 1 + math.exp(x[-1])

    return x, y, z, precise_ans


def runge_kutta_method(x_start, x_end, h):
    steps = int(round((x_end - x_start) / h))

    # приводим ур-ние к системе
    y = [0] * (steps + 1)
    z = [0] * (steps + 1)
    x = [x_start + i * h for i in range(steps + 1)]
    precise_ans = [0] * (steps + 1)

    y[0] = 2 + math.exp(1)
    z[0] = 1 + math.exp(1)

    for i in range(steps):
        precise_ans[i] = x[i] + 1 + math.exp(x[i])

        k_1 = h * func_y(y[i], z[i], h)
        l_1 = h * func_z(x[i], y[i], z[i], h)

        k_2 = h * func_y(y[i] + (k_1 / 2), z[i] + (l_1 / 2), h)
        l_2 = h * func_z(x[i] + (h / 2), y[i] + (k_1 / 2), z[i] + (l_1 / 2), h)

        k_3 = h * func_y(y[i] + (k_2 / 2), z[i] + (l_2 / 2), h)
        l_3 = h * func_z(x[i] + (h / 2), y[i] + (k_2 / 2), z[i] + (l_2 / 2), h)

        k_4 = h * func_y(y[i] + k_3, z[i] + l_3, h)
        l_4 = h * func_z(x[i] + h, y[i] + k_3, z[i] + l_3, h)

        y[i + 1] = y[i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        z[i + 1] = z[i] + (1 / 6) * (l_1 + 2 * l_2 + 2 * l_3 + l_4)

    precise_ans[-1] = x[-1] + 1 + math.exp(x[-1])

    return x, y, z, precise_ans


def adams_method(x_start, x_end, h):
    steps = int(round((x_end - x_start) / h))

    y = [0] * (steps + 1)
    z = [0] * (steps + 1)
    x = [x_start + i * h for i in range(steps + 1)]
    precise_ans = [0] * (steps + 1)

    # считаем первые 4 значения, чтобы составить интерполяционный полином
    x_r, y_r, z_r, precise_ans_r = runge_kutta_method(x_start, x_start + h * 3, h)
    for i in range(4):
        y[i] = y_r[i]
        z[i] = z_r[i]
        precise_ans[i] = precise_ans_r[i]

    for i in range(3, steps):
        y[i + 1] = y[i] + (h / 24) * (55 * func_y(y[i], z[i], h) - 59 * func_y(y[i - 1], z[i - 1], h) +
                                      37 * func_y(y[i - 2], z[i - 2], h) - 9 * func_y(y[i - 3], z[i - 3], h))

        z[i + 1] = z[i] + (h / 24) * (55 * func_z(x[i], y[i], z[i], h) - 59 * func_z(x[i- 1], y[i - 1], z[i - 1], h) +
                                      37 * func_z(x[i - 2], y[i - 2], z[i - 2], h) -
                                      9 * func_z(x[i - 3], y[i - 3], z[i - 3], h))
        precise_ans[i] = x[i] + 1 + math.exp(x[i])

    precise_ans[-1] = x[-1] + 1 + math.exp(x[-1])

    return x, y, z, precise_ans


def calculate_error(y_precise, y_real):
    errors = [0] * len(y_precise)
    for i in range(len(y_precise)):
        errors[i] = abs(y_precise[i] - y_real[i])
    return errors


def print_answer(x_start, x_end, h):
    h2 = h / 2
    _, y_2h, _, _ = euler_method(x_start, x_end, h2)

    x_e, y_e, z_e, precise = euler_method(x_start, x_end, h)
    errors_e = calculate_error(precise, y_e)

    print("============================================================================")
    print(f"метод Эйлера с шагом {h}:")
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    for i in range(len(x_e)):
        runge = runge_romberg_method(y_2h[2 * i], y_e[i], 1)
        print(f"{x_e[i]:<6.2f} | {y_e[i]:<15.6f} | {precise[i]:<15.6f} | {errors_e[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    _, y_2h, _, _ = runge_kutta_method(x_start, x_end, h2)

    x_r, y_r, z_r, _ = runge_kutta_method(x_start, x_end, h)
    errors_r = calculate_error(precise, y_r)

    print(f"\nметод Рунге-Кутты 4-го порядка с шагом {h}:")
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    for i in range(len(x_e)):
        runge = runge_romberg_method(y_2h[2 * i], y_r[i], 4)
        print(f"{x_r[i]:<6.2f} | {y_r[i]:<15.6f} | {precise[i]:<15.6f} | {errors_r[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    _, y_2h, _, _ = adams_method(x_start, x_end, h2)

    x_a, y_a, z_a, _ = adams_method(x_start, x_end, h)
    errors_a = calculate_error(precise, y_a)

    print(f"\nметод Адамса 4-го порядка с шагом {h}:")
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    for i in range(len(x_e)):
        runge = runge_romberg_method(y_2h[2 * i], y_a[i], 4)
        print(f"{x_a[i]:<6.2f} | {y_a[i]:<15.6f} | {precise[i]:<15.6f} | {errors_a[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    draw_graphics(x_start, x_end, h)


def draw_graphics(x_start, x_end, h):
    x_euler, y_euler, z_euler, exact_euler = euler_method(x_start, x_end, h)
    x_rk, y_rk, z_rk, exact_rk = runge_kutta_method(x_start, x_end, h)
    x_adams, y_adams, z_adams, exact_adams = adams_method(x_start, x_end, h)

    plt.figure(figsize=(10, 6))

    plt.scatter(x_euler, y_euler, color='blue', s=50, zorder=3, label='Метод Эйлера')
    plt.plot(x_euler, y_euler, color='blue', linewidth=2, zorder=2, label=f'Эйлер, h = {h}')

    plt.scatter(x_rk, y_rk, color='orange', s=50, zorder=3, label='Метод Рунге-Кутты 4-го порядка')
    plt.plot(x_rk, y_rk, color='orange', linewidth=2, zorder=2, alpha=0.5, label=f'Рунге-Кутта, h = {h}')

    plt.scatter(x_adams, y_adams, color='green', s=50, zorder=3, label='Метод Адамса 4-го порядка')
    plt.plot(x_adams, y_adams, color='green', linewidth=2, zorder=2, alpha=0.5, label=f'Адамс, h = {h}')

    plt.plot(x_rk, exact_rk, color='red', linewidth=2.5, zorder=1, label='Точное решение')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.show()


print_answer(1, 2, 0.1)
