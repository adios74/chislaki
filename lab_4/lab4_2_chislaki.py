import lab_1.lab2_chislaki as progonka
import matplotlib.pyplot as plt


def shooting_method(h, x_start, x_end):
    # начальные эта для задачи Коши
    eta_0 = -1.0
    eta_1 = -3.0

    epsilon = 0.0001

    # для начальных эта решаем задачу Коши методо Рунге-Кутты
    _, y_end_0, z_end_0 = runge_kutta_method(x_start, x_end, h, eta_0)
    phi_0 = phi(y_end_0[-1], z_end_0[-1])

    _, y_end_1, z_end_1 = runge_kutta_method(x_start, x_end, h, eta_1)
    phi_1 = phi(y_end_1[-1], z_end_1[-1])

    while abs(phi_1) > epsilon:
        # высчитывается новое приближение параметра методом секущих
        eta_next = eta_1 - ((eta_1 - eta_0) / (phi_1 - phi_0)) * phi_1

        # параметры сдвигаются
        eta_0 = eta_1
        phi_0 = phi_1
        eta_1 = eta_next

        # решается задача КОши для нового полученного значения эта
        _, y_end_1, z_end_1 = runge_kutta_method(x_start, x_end, h, eta_1)
        phi_1 = phi(y_end_1[-1], z_end_1[-1])

    # считаются все значения с полученным финальным параметром
    x, answers, _ = runge_kutta_method(x_start, x_end, h, eta_1)

    steps = int(round((x_end - x_start) / h))
    precise_ans = calculate_precise(steps, x)

    errors = [0.0] * (steps + 1)
    for i in range(steps + 1):
        errors[i] = abs(precise_ans[i] - answers[i])

    return x, answers, precise_ans, errors


def finite_difference_method(h, x_start, x_end):
    steps = int(round((x_end - x_start) / h))
    x = [x_start + i * h for i in range(steps + 1)]

    # считаем точное решение
    precise_ans = calculate_precise(steps, x)

    # инициализируем матрицу коэффициентов слау
    matrix = [[0.0] * (steps + 1) for _ in range(steps)]

    # левая граница
    p = (x[1] - 3) / (x[1] ** 2 - 1)
    q = -1 / (x[1] ** 2 - 1)
    f = 0
    a = 1 - (p * h) / 2
    b = -2 + h ** 2 * q
    c = 1 + (p * h) / 2
    matrix[0][0] = a + b
    matrix[0][1] = c

    # правая граница
    matrix[-1][steps - 2] = -1.0
    matrix[-1][steps - 1] = 1.0 + h
    matrix[-1][-1] = -0.75 * h

    # значения посередине
    for i in range(1, steps - 1):
        p = (x[i + 1] - 3) / (x[i + 1] ** 2 - 1)
        q = -1 / (x[i + 1] ** 2 - 1)
        f = 0
        matrix[i][i - 1] = 1 - (p * h) / 2
        matrix[i][i] = -2 + h ** 2 * q
        matrix[i][i + 1] = 1 + (p * h) / 2
        matrix[i][-1] = f * h ** 2

    # решаем полученную слау методом прогонки
    answers = calculate_system(matrix, steps)
    all_answers = [answers[0]] + answers

    # считаем ошибки
    errors = [0.0] * (steps + 1)
    for i in range(steps + 1):
        errors[i] = abs(precise_ans[i] - all_answers[i])

    return x, all_answers, precise_ans, errors


def runge_kutta_method(x_start, x_end, h, eta):
    steps = int(round((x_end - x_start) / h))

    y = [0] * (steps + 1)
    z = [0] * (steps + 1)
    h = (x_end - x_start) / steps
    x = [x_start + i * h for i in range(steps + 1)]

    y[0] = eta
    z[0] = 0.0

    for i in range(steps):
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

    return x, y, z


def func_y(y, z, h):
    return z


def func_z(x, y, z, h):
    if x >= 1.0 - 1e-9:
        return 0.25
    return (y - (x - 3) * z) / (x ** 2 - 1)


def phi(y_end, z_end):
    return y_end + z_end + 0.75


def calculate_precise(steps, x):
    precise_ans = [0] * (steps + 1)
    for i in range(steps + 1):
        precise_ans[i] = x[i] - 3 + 1 / (x[i] + 1)
    return precise_ans


def calculate_system(matrix, size):
    with open("test.txt", "w", encoding="utf-8") as file:
        file.write(str(size) + "\n")
        for i in range(size):
            for j in range(size + 1):
                file.write(str(matrix[i][j]) + " ")
            file.write("\n")

    matrix, size = progonka.read_elements("test.txt")
    answers, _ = progonka.tridiagonal_method(matrix, size)
    return answers


def runge_romberg_method(y_h, y_2h, p):
    return (y_h - y_2h) / (2 ** p - 1)


def print_answer(h1, h2, x_start, x_end):
    h_12 = h1 / 2
    _, y_2h, _, _ = shooting_method(h_12, x_start, x_end)
    print("============================================================================")
    print(f"метод стрельбы с шагом {h1}:")
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    x_s, answers_s, precise_ans_s, errors_s = shooting_method(h1, x_start, x_end)
    for i in range(len(x_s)):
        runge = runge_romberg_method(y_2h[2 * i], answers_s[i], 4)
        print(f"{x_s[i]:<6.2f} | {answers_s[i]:<15.6f} | {precise_ans_s[i]:<15.6f} | {errors_s[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    print(f"конечно-разностный метод с шагом {h1}:")
    _, y_2h, _, _ = finite_difference_method(h_12, x_start, x_end)
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    x_f, answers_f, precise_ans_f, errors_f = finite_difference_method(h1, x_start, x_end)
    for i in range(len(x_s)):
        runge = runge_romberg_method(y_2h[2 * i], answers_f[i], 2)
        print(f"{x_f[i]:<6.2f} | {answers_f[i]:<15.6f} | {precise_ans_f[i]:<15.6f} | {errors_f[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    print(f"метод стрельбы с шагом {h2}:")
    h_22 = h2 / 2
    _, y_2h, _, _ = shooting_method(h_22, x_start, x_end)
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    x_s, answers_s, precise_ans_s, errors_s = shooting_method(h2, x_start, x_end)
    for i in range(len(x_s)):
        runge = runge_romberg_method(y_2h[2 * i], answers_s[i], 4)
        print(f"{x_s[i]:<6.2f} | {answers_s[i]:<15.6f} | {precise_ans_s[i]:<15.6f} | {errors_s[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    print(f"конечно-разностный метод с шагом {h2}:")
    _, y_2h, _, _ = finite_difference_method(h_22, x_start, x_end)
    print(f"{'x':<6} | {'y':<15} | {'точное решение':<15} | {'ошибка':<15} | {'Рунге-Ромберг':<15}")
    print("----------------------------------------------------------------------------")
    x_f, answers_f, precise_ans_f, errors_f = finite_difference_method(h2, x_start, x_end)
    for i in range(len(x_s)):
        runge = runge_romberg_method(y_2h[2 * i], answers_f[i], 2)
        print(f"{x_f[i]:<6.2f} | {answers_f[i]:<15.6f} | {precise_ans_f[i]:<15.6f} | {errors_f[i]:<15.6e} | {runge:<15.6e}")
    print("============================================================================")

    draw_graphics(x_start, x_end, h1, h2)


def draw_graphics(x_start, x_end, h1, h2):
    x_s_1, y_s_1, _, _ = shooting_method(h1, x_start, x_end)
    x_f_1, y_f_1, _, _ = finite_difference_method(h1, x_start, x_end)

    x_s_2, y_s_2, precise, _ = shooting_method(h2, x_start, x_end)
    x_f_2, y_f_2, _, _ = finite_difference_method(h2, x_start, x_end)

    plt.figure(figsize=(10, 6))

    plt.scatter(x_s_1, y_s_1, color='blue', s=50, zorder=3, label=f'shooting method, step = {h1}')

    plt.scatter(x_f_1, y_f_1, color='orange', s=50, zorder=3, alpha=0.5, label=f'finite difference method, step = {h1}')

    plt.scatter(x_s_2, y_s_2, color='red', s=50, zorder=3, alpha=0.5, label=f'shooting method, step = {h2}')

    plt.scatter(x_f_2, y_f_2, color='green', s=50, zorder=3, alpha=0.5, label=f'finite difference method, step = {h2}')

    plt.plot(x_s_2, precise, color='black', linewidth=2.5, zorder=1, label='exact solution')

    plt.xlim(x_start - 0.05, x_end + 0.05)
    plt.ylim(-2.1, -1.3)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.show()


print_answer(0.1, 0.05, 0, 1)

