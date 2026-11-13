

def calculate_d(x, y, x_precise):
    y_left = (y[2] - y[1]) / (x[2] - x[1])
    y_right = (y[3] - y[2]) / (x[3] - x[2])

    y_1_d = y_left + ((y_right - y_left) / (x[3] - x[1])) * (2 * x_precise - x[1] - x[2])

    y_2_d = 2 * (y_right - y_left) / (x[3] - x[1])

    return y_left, y_right, y_1_d, y_2_d


def print_answer(x, y, x_precise, y_1_d, y_2_d):
    print(f"исходная таблица значений функций:")
    print(f"{'x':^10} | {'y':^10}")
    for i in range(len(x)):
        print(f"{x[i]:^10.2f} | {y[i]:^10.5f}")
    print(f"\nточка, в которой выполняются вычисления:\n{x_precise}\n")
    print(f"узлы таблицы, использованные при вычислении производной:")
    print(f"x[1] = {x[1]}, x[2] = {x[2]}, x[3] = {x[3]}\n")
    print(f"значение первой производной:\n{y_1_d}\n")
    print(f"значение второй производной:\n{y_2_d}\n")


x = [-0.2, 0.0, 0.2, 0.4, 0.6]
y = [-0.40136, 0.0, 0.40136, 0.81152, 1.2435]
x_precise = 0.2
y_left, y_right, y_1_d, y_2_d = calculate_d(x, y, x_precise)
print_answer(x, y, x_precise, y_1_d, y_2_d)