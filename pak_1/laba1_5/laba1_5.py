import math
import numpy as np


with open("input.txt", "r") as file_input:
    n = int(file_input.readline())
    eps = float(file_input.readline())

    A = []

    for _ in range(n):
        row = list(map(float, file_input.readline().split()))
        A.append(row)


    def qr(A):
        n = len(A)

        Q = [[0.0] * n for _ in range(n)]
        R = [[0.0] * n for _ in range(n)]

        for j in range(n): # столбец
            v = [0.0] * n
            for i in range(n):
                v[i] = A[i][j]

            for k in range(j): # r_12 = q_1^T * a2
                R[k][j] = 0.0

                for i in range(n):
                    R[k][j] += Q[i][k] * A[i][j]

                for i in range(n): # v = v - r12*q_1
                    v[i] -= R[k][j] * Q[i][k]

            norm = 0.0
            for i in range(n):
                norm += v[i] * v[i]

            norm = math.sqrt(norm)

            R[j][j] = norm

            if norm == 0:
                raise ValueError("Столбцы матрицы линейно зависимы")

            for i in range(n):
                Q[i][j] = v[i] / norm

        return Q, R


    def get_values(A, eps):
        n = len(A)
        values = []
        i = 0

        while i < n:
            # вещественное сз
            if i == n - 1 or abs(A[i + 1][i]) < eps:    # текущая строка последняя или элемент под диагональю почти 0
                values.append(A[i][i])
                i += 1

            # мб два св
            else:
                a = A[i][i]
                b = A[i][i + 1]
                c = A[i + 1][i]
                d = A[i + 1][i + 1]

                coef_b = a + d
                coef_c = a * d - b * c
                D = coef_b * coef_b - 4 * coef_c

                if D >= 0:
                    values.append((coef_b + math.sqrt(D)) / 2)
                    values.append((coef_b - math.sqrt(D)) / 2)
                else:
                    real = coef_b / 2
                    imag = math.sqrt(-D) / 2

                    values.append(complex(real, imag))
                    values.append(complex(real, -imag))

                i += 2

        return values

    A_original = [row[:] for row in A]
    Q_original, R_original = qr(A_original)
    previous_values = None

    while True:

        Q, R = qr(A)

        # A = R * Q
        new_A = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    new_A[i][j] += R[i][k] * Q[k][j]

        A = new_A

        values = get_values(A, eps)     # СЗ матрицы A_new
        has_complex = any(isinstance(value, complex) for value in values)

        if not has_complex:
            sum = 0.0

            for i in range(1, n):
                for j in range(i):
                    sum += A[i][j] ** 2

            if math.sqrt(sum) <= eps:
                break

        elif previous_values is not None:
            max_diff = max(
                abs(values[i] - previous_values[i])
                for i in range(len(values))
            )

            if max_diff <= eps:
                break
        previous_values = values[:]

    # для вывода в файл нашего разложения
    Q = Q_original
    R = R_original
    QR = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            for k in range(n):
                QR[i][j] += Q[i][k] * R[k][j]

    values = get_values(A, eps)

    # numpy
    B = np.array(A_original)
    complex_values = np.linalg.eigvals(B)

with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("1. Исходная матрица A:\n")
    for row in A_original:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n2. Матрица Q:\n")
    for row in Q:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n3. Матрица R:\n")
    for row in R:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n4. Результат Q * R:\n")
    for row in QR:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n5. Найденные собственные значения:\n")
    for value in values:
        file_output.write(f"{value} ")

    file_output.write("\n\n6. Матрица с комплексно-сопряжёнными собственными значениями:\n")
    for row in B:
        file_output.write(" ".join(f"{x}" for x in row) + "\n")

    file_output.write("\nСобственные значения:\n")
    for value in complex_values:
        file_output.write(f"{value}\n")
