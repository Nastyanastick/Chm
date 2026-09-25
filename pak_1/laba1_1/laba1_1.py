with open("input.txt", "r") as file_input:
    n = int(file_input.readline())
    A = []
    b = []
    det = 1.0
    swap_count = 0

    for _ in range(n):
        row = list(map(float, file_input.readline().split()))
        A.append(row[:n])
        b.append(row[n])

    L = [[0.0] * n for _ in range(n)]
    U = [A[i][:] for i in range(n)]
    P = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0
        P[i][i] = 1.0

    for k in range(n - 1):          # столбцы
        index_max = k               # индекс максимального элемента в стлб
        for i in range(k, n):       # строки
            if (abs(U[i][k]) > abs(U[index_max][k])):
                index_max = i
        if index_max != k:
            U[k], U[index_max] = U[index_max], U[k]         # поменяли местами.
            P[k], P[index_max] = P[index_max], P[k]
            for j in range(k):
                L[k][j], L[index_max][j] = L[index_max][j], L[k][j]
            swap_count += 1
        coeff = 1
        for i in range(k + 1, n):               # зануляем определенный столбец. пишем матрицу L.
            if (U[i][k] == U[k][k]):
                coeff = 1
            else:
                coeff = U[i][k] / U[k][k]
            L[i][k] = coeff
            for j in range(k, n):               # меняем всю строку, вычитаем.
                U[i][j] -= coeff * U[k][j]

    # перемножение LU.
    result = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i][j] += L[i][k] * U[k][j]

    # решение Ax=b
    y_1 = [0.0] * n
    pb = [0.0] * n

    for i in range(n):
        for j in range(n):
            pb[i] += P[i][j] * b[j]

    for i in range(n):
        sum = 0.0
        for j in range(i):
            sum += L[i][j] * y_1[j]
        y_1[i] = pb[i] - sum
        # y[i] = y[i] / L[i][i] нет смысла, тк диагональ L = 1

    # UX=y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum = 0.0
        for j in range(i + 1, n):
            sum += U[i][j] * x[j]
        x[i] = (y_1[i] - sum) / U[i][i]

    # решение Ax = I; (x_a = A^{-1})
    # LUx = I -> Ly = I(P)
    y_2 = [[0.0] * n for _ in range(n)]
    for k in range(n): # номер столбца y
        for i in range(n):
            sum = 0.0
            for j in range(i):
                sum += L[i][j] * y_2[j][k]
            y_2[i][k] = P[i][k] - sum

    # Ux=y
    x_a = [[0.0] * n for _ in range(n)]
    for k in range(n):
        for i in range(n - 1, -1, -1):
            sum = 0.0
            for j in range(i + 1, n):
                sum += U[i][j] * x_a[j][k]
            x_a[i][k] = (y_2[i][k] - sum) / U[i][i]

    # det(A) = det(L) det(U)
    det_u = 1.0
    for i in range(n):
        det_u *= U[i][i]
    if (swap_count % 2 == 0):
        det = det_u
    else:
        det = -det_u

    # проверка AA^{-1} = I
    check = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                check[i][j] += A[i][k] * x_a[k][j]

with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("1. L:\n")
    for row in L:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\nU:\n")
    for row in U:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n2. L * U:\n")
    for row in result:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n3. Ax = b:\n")
    file_output.write(" ".join(f"{x:.6f}" for x in x) + "\n")

    file_output.write("\n4. A^{-1}:\n")
    for row in x_a:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n5. det(A):\n")
    file_output.write(f"{det:.6f}\n")

    file_output.write("\n6. A * A^{-1} = E:\n")

    check_inverse = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            for k in range(n):
                check_inverse[i][j] += A[i][k] * x_a[k][j]

    for row in check_inverse:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\n7. L * U = P * A:\n")

    check_pa = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            for k in range(n):
                check_pa[i][j] += P[i][k] * A[k][j]

    file_output.write("L * U:\n")
    for row in result:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")

    file_output.write("\nP * A:\n")
    for row in check_pa:
        file_output.write(" ".join(f"{x:.6f}" for x in row) + "\n")