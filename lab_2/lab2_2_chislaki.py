import math

def simple_iteration_method(epsilon):
    x_1, x_2 = 0.25, 0.25

    d_f1_x = 0
    d_f1_y = - 1/3 * math.sin(x_2)
    d_f2_x = 1/3 * math.exp(x_1)
    d_f2_y = 0

    norm_1 = abs(d_f1_x) + abs(d_f1_y)
    norm_2 = abs(d_f2_x) + abs(d_f2_y)
    q = max(norm_1, norm_2)
    if q >= 1:
        print("не сходится")
        return False

    x_1_prev = 1 + x_1 - epsilon * (1 - q) / q
    x_2_prev = 0
    while abs((q / (1 - q)) * (x_1 - x_1_prev)) > epsilon or abs((q / (1 - q)) * (x_2 - x_2_prev)) > epsilon:
        x_1_prev = x_1
        x_2_prev = x_2
        x_1 = 1 / 3 * math.cos(x_2_prev)
        x_2 = 1 / 3 * math.exp(x_1_prev)

    return x_1, x_2


def newton_method(epsilon):
    x_1_prev, x_2_prev = 0, 0

    x_1, x_2 = 0.25, 0.25

    while abs(x_1 - x_1_prev) > epsilon or abs(x_2 - x_2_prev) > epsilon:
        f_1 = 3 * x_1 - math.cos(x_2)
        f_2 = 3 * x_2 - math.exp(x_1)

        d_f1_x = 3
        d_f1_y = math.sin(x_2)
        d_f2_x = - math.exp(x_1)
        d_f2_y = 3

        det_j = d_f1_x * d_f2_y - d_f2_x * d_f1_y
        det_a_1 = f_1 * d_f2_y - f_2 * d_f1_y
        det_a_2 = d_f1_x * f_2 - d_f2_x * f_1

        x_1_prev = x_1
        x_2_prev = x_2

        x_1 = x_1_prev - det_a_1 / det_j
        x_2 = x_2_prev - det_a_2 / det_j

    return x_1, x_2

print("\nметод Ньютона:")
print(newton_method(1e-4))
print("\nметод простых итераций:")
print(simple_iteration_method(1e-4))