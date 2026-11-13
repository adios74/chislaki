import math
import numpy as np

def read_matrix(filename):
    try:
        with open(filename, "r", encoding="utf-8") as file:
            epsilon = float(file.readline().strip())
            size = int(file.readline().strip())
            matrix = [[0] * size for _ in range(size)]
            for i in range(size):
                numbers = file.readline().split()
                for j in range(size):
                    matrix[i][j] = float(numbers[j])

    except FileNotFoundError:
        print("файл не найден")
        return None, None, None
    return size, epsilon, matrix


def phi_calculation(matrix, line, column):
    if matrix[line][line] == matrix[column][column]:
        return math.pi / 4

    return 0.5 * math.atan(2 * matrix[line][column] / (matrix[line][line] - matrix[column][column]))


def rotation_matrix(size, line, column, sin_phi, cos_phi):
    U = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if i == line and j == column:
                U[i][j] = -sin_phi
            elif i == column and j == line:
                U[i][j] = sin_phi
            elif (i == column and j == column) or (i == line and j == line):
                U[i][j] = cos_phi
            elif (i == j):
                U[i][j] = 1
            else:
                U[i][j] = 0
    return U


def check_result(diagonal, matrix, size, eigenvectors, epsilon):
    A = np.array(matrix)
    V = np.array(eigenvectors)

    Lambda = [[0] * size for _ in range(size)]
    counter = 0
    for i in range(size):
        for j in range(size):
            if i == j:
                Lambda[i][j] = diagonal[counter]
                counter += 1
    Lambda = np.array(Lambda)
    if np.allclose(A @ V, V @ Lambda, atol=epsilon):
        return True
    else:
        return False


def print_result(matrix, size, epsilon):
    diagonal, eigenvectors = yakobi(matrix, size, epsilon)

    if not check_result(diagonal, matrix, size, eigenvectors, epsilon):
        print("не пройдена проверка\n")
        return False

    print("проверка пройдена!\n")
    print(f"заданная точность вычислений:\n{epsilon}\n")

    print("найденные СЗ:")
    for i in range(size):
        print(f"СЗ{i} = {diagonal[i]}", end="  ")

    print("\n\nнайденные СВ:")
    counter = 0
    while counter < size:
        print("[", end=" ")
        for i in range(size):
            for j in range(size):
                if j == counter:
                    print(f" {eigenvectors[i][j]} ", end=" ")
        print("] ^ T\n")
        counter += 1

    print("найденная матрица СВ:")
    for i in range(size):
        print("\n", end=" ")
        for j in range(size):
            print(f" {eigenvectors[i][j]} ", end=" ")
    print()
    return True


def yakobi(matrix, size, epsilon):

    t = epsilon + 1
    diagonal = [0] * size

    # составление матрицы собственных векторов
    eigenvectors = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if i == j:
                eigenvectors[i][j] = 1
    eigenvectors = np.array(eigenvectors)

    while t > epsilon:
        # 1. поиск наибольшего по модулю внедиагонального элемента
        maximum = 0
        line = 0
        column = 0
        for i in range(size):
            for j in range(size):
                if i != j:
                    if abs(matrix[i][j]) > maximum:
                        maximum = abs(matrix[i][j])
                        line = i
                        column = j

        # 2. нахождение матрицы вращения
        phi = phi_calculation(matrix, line, column)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)

        U = rotation_matrix(size, line, column, sin_phi, cos_phi)
        eigenvectors = eigenvectors @ U
        transpositioned_U = np.transpose(U)
        U = np.array(U)
        transpositioned_U = np.array(transpositioned_U)

        # 3. перемножение полученных матриц с исходной
        matrix = transpositioned_U @ matrix @ U

        # 4. проверка вычисленных значений
        sum_not_diagonal = 0
        for l in range(size):
            for m in range(size):
                if l < m:
                    sum_not_diagonal += matrix[l][m] ** 2
                if l == m:
                    diagonal[l] = matrix[l][m]
        t = math.sqrt(sum_not_diagonal)

    return diagonal, eigenvectors


print("введите название файла для теста:")
filename = input()
size, epsilon, matrix = read_matrix(filename)
print_result(matrix, size, epsilon)