with open("input.txt") as file_input:
    n = int(file_input.readline())
    eps = float(file_input.readline())
    A = []
    b = []
    iterators_simple = 0
    iterators_seidel = 0

    # система Ax=b
    for _ in range(n):
        row = list(map(float, file_input.readline().split()))

        A.append(row[:-1])
        b.append(row[-1])

    x = [0.0] * n       # начальное приближение

    while True:
        new_x = [0.0] * n

        for i in range(n):
            sum = 0.0
            for j in range(n):
                if j != i:
                    sum += A[i][j] * x[j]
            new_x[i] = (b[i] - sum) / A[i][i]

        # точность
        error = 0.0

        for i in range(n):
            if abs(new_x[i] - x[i]) > error:
                error = abs(new_x[i] - x[i])
        x = new_x
        iterators_simple += 1

        if error < eps:
            break

    x_simple = x[:]
    x = [0.0] * n

    while True:
        old_x = x[:]        # сохраненное старое

        for i in range(n):
            sum = 0.0

            for j in range(n):
                if j != i:
                    sum += A[i][j] * x[j]
            x[i] = (b[i] - sum) / A[i][i]

        error = 0.0

        for i in range(n):
            if abs(x[i] - old_x[i]) > error:
                error = abs(x[i] - old_x[i])

        iterators_seidel += 1

        if error < eps:
            break

    x_seidel = x[:]

    # метод гаусса для эталонного решения

    A_gauss = [A[i][:] for i in range(n)]
    b_gauss = b[:]

    # прямой ход
    for k in range(n - 1):
        index_max = k
        for i in range(k + 1, n):
            if abs(A_gauss[i][k]) > abs(A_gauss[index_max][k]):
                index_max = i
        if index_max != k:
            A_gauss[k], A_gauss[index_max] = A_gauss[index_max], A_gauss[k]
            b_gauss[k], b_gauss[index_max] = b_gauss[index_max], b_gauss[k]
        for i in range(k + 1, n):
            coeff = A_gauss[i][k] / A_gauss[k][k]
            for j in range(k, n):
                A_gauss[i][j] -= coeff * A_gauss[k][j]
            b_gauss[i] -= coeff * b_gauss[k]

    # обратный ход
    x_gauss = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum = 0.0
        for j in range(i + 1, n):
            sum += A_gauss[i][j] * x_gauss[j]
        x_gauss[i] = (b_gauss[i] - sum) / A_gauss[i][i]

    with open("output.txt", "w", encoding="utf-8") as file_output:

        file_output.write("Заданная точность вычислений:\n")
        file_output.write(f"{eps}\n")

        file_output.write("\nМетод простых итераций:\n")
        file_output.write("Решение:\n")
        for i in range(n):
            file_output.write(f"x{i + 1} = {x_simple[i]:.6f}\n")
        file_output.write(f"Количество итераций: {iterators_simple}\n")

        file_output.write("\nМетод Зейделя:\n")
        file_output.write("Решение:\n")
        for i in range(n):
            file_output.write(f"x{i + 1} = {x_seidel[i]:.6f}\n")
        file_output.write(f"Количество итераций: {iterators_seidel}\n")

        file_output.write("\nЭталонное решение:\n")
        for i in range(n):
            file_output.write(f"x{i + 1} = {x_gauss[i]:.6f}\n")

        file_output.write("\nСравнение методов:\n")
        if iterators_simple < iterators_seidel:
            file_output.write("Метод простых итераций сошелся за меньшее количество итераций.\n")
        elif iterators_seidel < iterators_simple:
            file_output.write("Метод Зейделя сошелся за меньшее количество итераций.\n")
        else:
            file_output.write("Оба метода сошлись за одинаковое количество итераций.\n")