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
    return epsilon, size, matrix


def qr_method(matrix, size):

    E = [[0] * size for _ in range(size)]
    Q = [[0] * size for _ in range(size)]
    for i in range(size):
        E[i][i] = 1
        Q[i][i] = 1
    E = np.array(E, dtype=float)
    Q = np.array(Q, dtype=float)
    R = np.array(matrix, dtype=float)

    for iteration in range(size - 1):
        # создание массива столбцов
        matrix_columns = [[0] * size for _ in range(size)]
        for i in range(size):
            for j in range(size):
                matrix_columns[j][i] = R[i][j]

        # 1. поиск вектора v
        norm = 0
        for k in range(iteration, size):
            norm += matrix_columns[iteration][k] ** 2
        norm = math.sqrt(norm)
        if norm == 0:
            continue

        vector = [0] * size
        for i in range(size):
            if i < iteration:
                vector[i] = 0
            elif i == iteration:
                if R[i][i] == 0:
                    vector[i] = R[i][i] + norm
                else:
                    vector[i] = R[i][i] + np.sign(R[i][i]) * norm
            elif i > iteration:
                vector[i] = R[i][iteration]

        # 2. составление матрицы Хаусхолдера
        v = np.array(vector)
        v_column = v.reshape(-1, 1)
        v_row = v.reshape(1, -1)
        v_chisl = v_column @ v_row
        v_znam = v @ v
        H = E - 2 * v_chisl / v_znam
        H = np.array(H)

        Q = Q @ H
        R = H @ R

    return Q, R


def calculate_eigenvalues(A, size, epsilon, maximum_iterations=100):
    # 3. итерационный процесс
    A = np.array(A, dtype=float)
    for count in range(maximum_iterations):
        ok = True
        Q, R = qr_method(A, size)
        A = R @ Q

        for i in range(1, size):
            for j in range(i):
                if abs(A[i][j]) >= abs(epsilon):
                    ok = False
                    break
        if ok:
            break

    # 4. нахождение собственных значений
    eigenvalues = []
    i = 0
    while i < size:
        if i == size - 1 or abs(A[i + 1][i]) < abs(epsilon):
            eigenvalues.append(A[i][i])
            i += 1
        else:
            a = A[i][i]
            b = A[i][i + 1]
            c = A[i + 1][i]
            d = A[i + 1][i + 1]
            trace = a + d
            det = a * d - b * c
            discriminant = trace ** 2 - 4 * det
            if discriminant >= 0:
                sq_root = discriminant ** 0.5
                eigenvalues.append((trace - sq_root) / 2)
                eigenvalues.append((trace + sq_root) / 2)
            else:
                discriminant *= (-1)
                sq_root = discriminant ** 0.5
                eigenvalues.append(complex(trace / 2, sq_root / 2))
                eigenvalues.append(complex(trace / 2, -sq_root / 2))
            i += 2

    return eigenvalues


def print_result(matrix, size, epsilon):
    A = np.array(matrix, dtype=float)
    print(f"исходная матрица A:\n{A}\n")

    Q, R = qr_method(A, size)
    print(f"матрица Q:\n{Q}\n")
    print(f"матрица R:\n{R}\n")

    multiplication = Q @ R
    print(f"произведение Q и R:\n{multiplication}\n")

    print("фактические собственные значения:")
    eigenvalues = calculate_eigenvalues(A, size, epsilon)
    for i in eigenvalues:
        print(i)

    print("\nэталонные собственные значения:")
    reference = np.linalg.eigvals(A)
    for i in reference:
        print(i)


print("введите название файла для теста:")
filename = input()
epsilon, size, matrix = read_matrix(filename)
print_result(matrix, size, epsilon)