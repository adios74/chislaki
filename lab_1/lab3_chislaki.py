import copy
import numpy as np


def simple_iteration(matrix, epsilon, size):
    matrix = copy.deepcopy(matrix)

    # 1. посчет первых приближений
    approximation = []
    for k in range(size):
        if matrix[k][k] == 0:
            approximation.append(0)
            continue
        approximation.append(matrix[k][-1] / matrix[k][k])

    # 2. выражение неизвестной из уравнения
    for i in range(size):
        for j in range(size):
            if i == j:
                continue
            matrix[i][j] /= matrix[i][i]
            matrix[i][j] *= -1
        matrix[i][size] /= matrix[i][i]

    # 3. подсчет значений с подстановкой предыдущих приближений и проверка условий
    iteration_count = 0
    ok = False
    while ok is False:
        iteration_count += 1
        ok = True
        current = []
        for i in range(size):
            x = 0
            for j in range(size):
                if i != j:
                    x += matrix[i][j] * approximation[j]
            x += matrix[i][size]
            current.append(x)
        for k in range(size):
            delta = (current[k] - approximation[k]) / current[k]
            if abs(delta) >= abs(epsilon):
                ok = False
        approximation = current.copy()

    return approximation, iteration_count


def zeidel_method(matrix, epsilon, size):
    matrix = copy.deepcopy(matrix)

    # 1. посчет первых приближений
    approximation = []
    for k in range(size):
        if matrix[k][k] == 0:
            approximation.append(0)
            continue
        approximation.append(matrix[k][-1] / matrix[k][k])

    # 2. выражение неизвестной из уравнения
    for i in range(size):
        for j in range(size):
            if (i == j):
                continue
            matrix[i][j] /= matrix[i][i]
            matrix[i][j] *= -1
        matrix[i][size] /= matrix[i][i]

    # 3. подсчет значений с подстановкой предыдущих приближений и проверка условий
    iteration_count = 0
    ok = False
    while ok is False:
        iteration_count += 1
        ok = True
        current = []
        for i in range(size):
            x = 0
            for j in range(size):
                if j > i:
                    x += matrix[i][j] * approximation[j]
                elif j < i:
                    x += matrix[i][j] * current[j]
            x += matrix[i][size]
            current.append(x)
        for k in range(size):
            delta = (current[k] - approximation[k]) / current[k]
            if abs(delta) >= abs(epsilon):
                ok = False
        approximation = current.copy()

    return approximation, iteration_count


def read_matrix(filename):
    try:
        with open(filename, "r", encoding="utf-8") as file:
            epsilon = float(file.readline().strip())
            size = int(file.readline().strip())
            matrix = [[0] * (size + 1) for _ in range(size)]
            for i in range(size):
                numbers = file.readline().split()
                for j in range(size + 1):
                    matrix[i][j] = float(numbers[j])

    except FileNotFoundError:
        print("файл не найден")
        return None, None, None
    return size, epsilon, matrix


def print_answer(matrix, epsilon, size):

    # 1. Заданную точность вычислений ε.
    # 2. Количество итераций для каждого метода.
    # 3. Точное или эталонное решение для проверки результата.
    # 4. Сравнение методов по скорости сходимости.

    approximation_s, iteration_count_s = simple_iteration(matrix, epsilon, size)
    approximation_z, iteration_count_z = zeidel_method(matrix, epsilon, size)

    print(f"заданная точность вычислений:\n{epsilon}\n")
    print(f"количество итераций для метода простых итераций:\n{iteration_count_s}\n")
    print(f"количество итераций для метода Зейделя:\n{iteration_count_z}\n")

    print(f"результат для метода простых итераций:\n{approximation_s}\n")
    print(f"результат для метода Зейделя:\n{approximation_z}\n")

    matrix_for_ref = [[0] * size for _ in range(size)]
    answers_for_ref = []
    for i in range(size):
        for j in range(size + 1):
            if j != size:
                matrix_for_ref[i][j] = matrix[i][j]
            else:
                answers_for_ref.append(matrix[i][j])

    matrix_for_ref = np.array(matrix_for_ref)
    reference = np.linalg.solve(matrix_for_ref, answers_for_ref)
    print(f"эталонное решение:\n{reference}\n")

    if iteration_count_s < iteration_count_z:
        print("метод простых итераций на данном примере сходится быстрее")
    elif iteration_count_s > iteration_count_z:
        print("метод Зейделя на данном примере сходится быстрее")
    else:
        print("на данном примере методы сходятся одинаково быстро")

    return True


# print("введите название файла для теста:")
# filename = input()
# size, epsilon, matrix = read_matrix(filename)
# print_answer(matrix, epsilon, size)

