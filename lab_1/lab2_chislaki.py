
def read_elements(filename):
    try:
        with open(filename, "r", encoding="utf-8") as file:
            size = int(file.readline().strip())
            elements = []
            for i in range(size):
                numbers = [float(num) for num in file.readline().split()]
                answers = numbers[-1]
                if i == 0:
                    raw = [numbers[0], numbers[1], answers]
                elif i == size - 1:
                    raw = [numbers[size - 2], numbers[size - 1], answers]
                else:
                    raw = [numbers[i - 1], numbers[i], numbers[i + 1], answers]
                elements.append(raw)

    except FileNotFoundError:
        print("файл не найден")
        return None, None
    return elements, size


def tridiagonal_method(elements, size):
    coefficients = [[0] for i in range(size)]
    element_num = 0
    last_x = 0
    answers = [0] * size

    for i in range(size):
        # преобразуем элементы строки, выражаем нужный элемент
        iterations = len(elements[i])
        for j in range(iterations):
            if j != element_num:
                elements[i][j] /= elements[i][element_num]
            if j != iterations - 1 and j != element_num:
                elements[i][j] *= -1
        elements[i][element_num] = 1

        if i == size - 1:
            last_x = elements[i][-1]
            coefficients[i] = [0, last_x]
            break

        # получаем прогоночные коэффициенты
        P = elements[i][element_num + 1]
        Q = elements[i][element_num + 2]
        coefficients[i] = [P, Q]

        # в следующем уравнении подставляем полученные прогоночные коэффициенты
        elements[i + 1][1] += elements[i + 1][0] * P
        elements[i + 1][-1] -= elements[i + 1][0] * Q
        elements[i + 1][0] = 0

        element_num = 1

    answers[size - 1] = last_x
    for i in range(size - 2, -1, -1):
        answers[i] = coefficients[i][0] * last_x + coefficients[i][1]
        last_x = answers[i]

    return answers, coefficients


def print_answer(elements, size):
    answers, coefficients = tridiagonal_method(elements, size)
    for i in range(size):
        print(f"x{i + 1} = {answers[i]}, P = {coefficients[i][0]}, Q = {coefficients[i][1]}")


# print("введите название файла для теста:")
# filename = input()
# elements, size = read_elements(filename)
# print_answer(elements, size)