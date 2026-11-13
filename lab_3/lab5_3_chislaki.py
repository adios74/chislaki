def calculate_func(x):
    return 1 / (256 - x ** 4)


def rectangle_method(x_0, x_k, h):
    steps = int(round((x_k - x_0) / h))

    f_left, f_right, f_centre = 0, 0, 0
    x_p, x_n = x_0, x_0 + h

    for i in range(steps):
        # метод средних прямоугольников
        x_c = (x_p + x_n) / 2
        f_centre += calculate_func(x_c)

        # метод левых прямоугольников
        f_left += calculate_func(x_p)

        # метод правых прямоугольников
        f_right += calculate_func(x_n)
        x_p = x_n
        x_n += h

    f_left *= h
    f_right *= h
    f_centre *= h

    return f_left, f_right, f_centre


def trapezoid_method(x_0, x_k, h):
    steps = int(round((x_k - x_0) / h))
    # крайние точки
    f = calculate_func(x_0) / 2 + calculate_func(x_k) / 2

    # подсчет площадей трапеций между точками
    x_cur = x_0 + h
    for i in range(1, steps):
        f += calculate_func(x_cur)
        x_cur += h

    f *= h
    return f


def simpson_method(x_0, x_k, h):
    steps = int(round((x_k - x_0) / h))

    if steps % 2 != 0:
        print("число шагов должно быть четным!")

    # крайние точки
    f = calculate_func(x_0) + calculate_func(x_k)

    x = x_0 + h
    for i in range(1, steps):
        if i % 2 == 0:
            f += 2 * calculate_func(x)
        else:
            f += 4 * calculate_func(x)
        x += h

    f *= h / 3
    return f


def runge_romberg_method(f_h, f_k_h, k, p):
    return f_h + (f_h - f_k_h) / (k ** p - 1)


def print_answer(x_0, x_k, h_1, h_2, precise_ans):
    k = int(h_1 / h_2)
    f_left_1, f_right_1, f_centre_1 = rectangle_method(x_0, x_k, h_1)
    f_1_t = trapezoid_method(x_0, x_k, h_1)
    f_1_s = simpson_method(x_0, x_k, h_1)

    f_left_2, f_right_2, f_centre_2 = rectangle_method(x_0, x_k, h_2)
    f_2_t = trapezoid_method(x_0, x_k, h_2)
    f_2_s = simpson_method(x_0, x_k, h_2)

    methods_data = [
        ("левые прямоугольники", f_left_1, f_left_2, 1),
        ("правые прямоугольники", f_right_1, f_right_2, 1),
        ("средние прямоугольники", f_centre_1, f_centre_2, 2),
        ("трапеции", f_1_t, f_2_t, 2),
        ("Симпсон", f_1_s, f_2_s, 4)
    ]

    print(f"{'метод':<24} | {'шаг h1':^10} | {'шаг h2':^10} | {'порядок точности':^10} | {'метод Рунге-Ромберга':^20} | {'погрешность':^12}")

    runge_all = []
    for method, f_1, f_2, p in methods_data:
        runge = runge_romberg_method(f_2, f_1, k, p)
        runge_all.append(runge)
        error = abs(runge - f_2)
        print(f"{method:<24} | {f_1:>10.5f} | {f_2:>10.5f} | {p:^16} | {runge:>20.5f} | {error:>12.2e}")

    print(f"аналитически полученное решение:\n{precise_ans}\n")
    errors = []
    errors.append(abs(runge_all[0] - precise_ans))
    errors.append(abs(runge_all[1] - precise_ans))
    errors.append(abs(runge_all[2] - precise_ans))
    errors.append(abs(runge_all[3] - precise_ans))
    errors.append(abs(runge_all[4] - precise_ans))
    print(f"абсолютная погрешность каждого метода:\n")
    print(f"левые прямоугольники - {errors[0]}, правые прямоугольники - {errors[1]}, средние прямоугольники - {errors[2]},"
          f"трапеции - {errors[3]}, Симпсон - {errors[4]}")

    min_err = float("inf")
    min_err_index = 6
    for i in range(5):
        if errors[i] < min_err:
            min_err = errors[i]
            min_err_index = i

    print(f"метод с наиболее точным результатом - {methods_data[min_err_index][0]}")


print_answer(-2.0, 2.0, 1.0, 0.5, 0.015828)