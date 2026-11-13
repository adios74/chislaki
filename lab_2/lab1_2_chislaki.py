import math


def newton_method(epsilon):
    x_prev = 2
    approximation = epsilon + 1
    x = 0

    while abs(approximation) > abs(epsilon):
        f_prev = 4 ** x_prev - 5 * x_prev - 2
        derivative_prev = 4 ** x_prev * math.log(4) - 5
        x = x_prev - (f_prev / derivative_prev)
        approximation = x_prev - x
        x_prev = x

    return x


def simple_iteration_method(epsilon):
    approximation = epsilon + 1
    x_prev = 1.5
    x = 0

    derivative_1 = (4 ** x_prev * math.log(4)) / 5
    derivative_2 = 5 / ((5 * x_prev + 2) * math.log(4))

    if abs(derivative_1) < 1:
        while abs(approximation) > abs(epsilon):
            x = (4 ** x_prev - 2) / 5
            approximation = x - x_prev
            x_prev = x

    elif abs(derivative_2) < 1:
        while abs(approximation) > abs(epsilon):
            x = math.log((5 * x_prev + 2), 4)
            approximation = x - x_prev
            x_prev = x

    return x

print("\nметод Ньютона:")
print(newton_method(1e-4))
print("\nметод простых итераций:")
print(simple_iteration_method(1e-4))