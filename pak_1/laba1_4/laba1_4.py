import math


with open("input.txt") as file_input:
    n = int(file_input.readline())
    eps = float(file_input.readline())

    A = []

    for _ in range(n):
        row = list(map(float, file_input.readline().split()))
        A.append(row)

    Q = [[0.0] * n for _ in range(n)]

    for i in range(n):
        Q[i][i] = 1.0

    A_original = [row[:] for row in A]

    while True:
        sum = 0.0
        max_value = 0.0
        p = 0
        q = 1

        for i in range(n):
            for j in range(i + 1, n):
                sum += A[i][j] ** 2
                if abs(A[i][j]) > max_value:
                    max_value = abs(A[i][j])
                    p = i
                    q = j

        if math.sqrt(sum) < eps:
            break

        # угол вращения
        phi = 0.5 * math.atan2(
            2.0 * A[p][q],
            A[p][p] - A[q][q]
        )

        c = math.cos(phi)
        s = math.sin(phi)

        R = [[0.0] * n for _ in range(n)]

        for i in range(n):
            R[i][i] = 1.0

        R[p][p] = c
        R[q][q] = c
        R[p][q] = -s
        R[q][p] = s

        RT = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                RT[i][j] = R[j][i]

        # R^T * A
        RTA = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    RTA[i][j] += RT[i][k] * A[k][j]

        # (R^T * A) * R
        new_A = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    new_A[i][j] += RTA[i][k] * R[k][j]
        A = new_A

        # Q *R
        QR = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    QR[i][j] += Q[i][k] * R[k][j]

        Q = QR

    Values = [0.0] * n
    for i in range(n):
        Values[i] = A[i][i]

    # проверка A_original * Q = Q * Matrix_values
    AQ = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            for k in range(n):
                AQ[i][j] += A_original[i][k] * Q[k][j]

    QValues = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            QValues[i][j] = Q[i][j] * Values[j]

with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("1. Заданная точность вычислений:\n")
    file_output.write(f"{eps:.6f}\n")

    file_output.write("\n2. Найденные собственные значения:\n")
    file_output.write(" ".join(f"{x:.6f}" for x in Values) + "\n")

    file_output.write("\n3. Найденные собственные векторы:\n")
    for i in range(n):
        file_output.write(
            " ".join(f"{Q[j][i]:.6f}" for j in range(n)) + "\n"
        )

    file_output.write("\nМатрица собственных векторов Q:\n")
    for row in Q:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n4. Проверка A * Q = Q * Λ:\n")

    file_output.write("A * Q:\n")
    for row in AQ:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\nQ * Λ:\n")
    for row in QValues:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")