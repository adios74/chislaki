import numpy as np


def find_main_element(matrix, size, i, P, b, permutation_counter):
    matrix = np.array(matrix)
    # поиск максимального ведущего элемента
    maximum = i
    for m in range(i + 1, size):
        if abs(matrix[m][i]) > abs(matrix[maximum][i]):
            maximum = m

    if maximum == i:
        return matrix[i][i], matrix, P, maximum, b, permutation_counter

    # составление матрицы для перестановки строк
    permutation = [[0] * size for i in range(size)]
    for m in range(size):
        for k in range(size):
            if m == i and k == maximum:
                permutation[m][k] = 1
            elif m == maximum and k == i:
                permutation[m][k] = 1
            elif m == k and m != i and m != maximum:
                permutation[m][k] = 1
    permutation_counter += 1
    permutation = np.array(permutation)
    P = permutation @ P
    matrix = permutation @ matrix
    b = permutation @ b
    return matrix[i][i], matrix, P, maximum, b, permutation_counter


def lu_method_forward(matrix, size, b):
    b = np.array(b, dtype=float)
    permutation_counter = 0

    # составление начальной матрицы перестановок (единичная матрица)
    P = [[0] * size for i in range(size)]
    for i in range(size):
        for j in range(size):
            if i == j:
                P[i][j] = 1
    P = np.array(P)

    L = [[0] * size for i in range(size)]
    for i in range(size):
        for j in range(size):
            if i == j:
                L[i][j] = 1

    U = np.array(matrix, dtype=float).copy()

    iteration = 0
    while iteration != size:
        # получение начального ведущего элемента
        main_element, U, P, maximum, b, permutation_counter = find_main_element(U, size, iteration, P, b, permutation_counter)

        # при изменении порядка строк меняем его и в матрице L
        if maximum != iteration:
            for j in range(iteration):
                L[iteration][j], L[maximum][j] = L[maximum][j], L[iteration][j]

        # проход по строкам матрицы и приведение к ступенчатому виду
        for i in range(iteration + 1, size):
            division_coef = U[i][iteration] / main_element
            L[i][iteration] = division_coef
            b[i] -= division_coef * b[iteration]
            for j in range(size):
                U[i][j] -= division_coef * U[iteration][j]
        iteration += 1

    return P, L, U, b, permutation_counter


def lu_method_backward_vector(size, U, b):
    answers = [0.0] * size
    for i in range(size - 1, -1, -1):
        x = 0
        for j in range(size - 1, -1, -1):
            if j != i:
                x -= answers[j] * U[i][j]
        x += b[i]
        x /= U[i][i]
        answers[i] = x
    return answers


def lu_method_backward_matrix(size, U, b):
    answers = [[0] * size for i in range(size)]
    answers = np.array(answers, dtype=float)
    for k in range(size):
        answers[:, k] = lu_method_backward_vector(size, U, b[:, k])
    return answers


def calculate_det(U, size, permutation_counter):
    det = 1
    for i in range(size):
        det *= U[i][i]
    det *= (-1) ** permutation_counter
    return det


def calculate_reversed(matrix, size):
    B = [[0] * size for i in range(size)]
    for i in range(size):
        for j in range(size):
            if i == j:
                B[i][j] = 1

    P, L, U, B, permutation_counter = lu_method_forward(matrix, size, B)
    reversed_matrix = lu_method_backward_matrix(size, U, B)
    return reversed_matrix


def check_result(A, A_reversed, L, U, P, b, size):
    epsilon = 1e-6
    A = np.array(A, dtype=float)
    A_reversed = np.array(A_reversed, dtype=float)
    reverse_check = A @ A_reversed

    reverse_ok = True
    for i in range(size):
        for j in range(size):
            if i == j:
                if abs(reverse_check[i][j] - 1.0) > epsilon:
                    reverse_ok = False
            elif i != j:
                if abs(reverse_check[i][j]) > epsilon:
                    reverse_ok = False

    lupa_ok = False
    L = np.array(L, dtype=float)
    U = np.array(U, dtype=float)
    P = np.array(P, dtype=float)

    if np.allclose(L @ U, P @ A, atol=epsilon):
        lupa_ok = True

    return lupa_ok, reverse_ok


def print_result(A, size, b):
    P, L, U, b, permutation_counter = lu_method_forward(A, size, b)

    L = np.array(L, dtype=float)
    U = np.array(U, dtype=float)
    print(f"матрица L:\n{L}\n")
    print(f"матрица U:\n{U}\n")

    multiplication = L @ U
    print(f"произведение матриц L и U:\n{multiplication}\n")

    answers = lu_method_backward_vector(size, U, b)
    print("решения системы Ax = b:")
    for i in range(len(answers)):
        print(f"x{i} = {answers[i]}")

    reversed_matrix = calculate_reversed(A, size)
    print(f"\nобратная матрица:\n{reversed_matrix}\n")

    det = calculate_det(U, size, permutation_counter)
    print(f"определитель матрицы:\n{det}\n")

    lupa_ok, reverse_ok = check_result(A, reversed_matrix, L, U, P, b, size)

    if not lupa_ok:
        print("ОШИБКА! выражение L * U = P * A не прошло проверку!")
        return False
    if not reverse_ok:
        print("ОШИБКА! выражение A * A ^ (-1) не прошло проверку!")
        return False
    print("все проверки пройдены!\nL * U = P * A\tи\tA * A ^ (-1)")
    return True


def read_matrix(filename):
    try:
        with open(filename, "r", encoding="utf-8") as file:
            size = int(file.readline().strip())
            matrix = [[0.0] * size for _ in range(size)]
            b = [0.0] * size
            for i in range(size):
                numbers = file.readline().split()
                for j in range(size):
                    matrix[i][j] = float(numbers[j])
                b[i] = float(numbers[size])

    except FileNotFoundError:
        print("файл не найден")
        return None, None, None
    return size, matrix, b


# print("введите название файла для теста:")
# filename = input()
# size, matrix, b = read_matrix(filename)
# print_result(matrix, size, b)