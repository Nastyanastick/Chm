with open("input.txt") as file_input:
    n = int(file_input.readline())

    a = [0.0] * (n - 1)     # нижняя
    b = [0.0] * n       # главная
    c = [0.0] * (n - 1)      # верхняя
    d = [0.0] * n
    x = [0.0] * n

    for i in range(n):
        row = list(map(float, file_input.readline().split()))

        if i == 0:
            b[i] = row[0]
            c[i] = row[1]
            d[i] = row[2]
        elif i == n - 1:
            a[i - 1] = row[0]
            b[i] = row[1]
            d[i] = row[2]
        else:
            a[i - 1] = row[0]
            b[i] = row[1]
            c[i] = row[2]
            d[i] = row[3]

    alpha = [0.0] * n
    beta = [0.0] * n
    alpha[0] = d[0] / b[0]
    beta[0] = c[0] / b[0]

    for i in range(1, n):
        denomirator = b[i] - a[i - 1] * beta[i - 1]
        alpha[i] = (d[i] - a[i - 1] * alpha[i - 1]) / denomirator
        if i < n - 1:       # последняя строка нам уже вся известна
            beta[i] = c[i] / denomirator

    x[n - 1] = alpha[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] - beta[i] * x[i + 1]

with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("Коэффициенты alpha: ")
    file_output.write(" ".join(f"{x:.6f}" for x in alpha) + "\n")
    file_output.write("Коэфициенты beta: ")
    file_output.write(" ".join(f"{x:.6f}" for x in beta) + "\n")

    file_output.write("Решение системы уравнений:\n")
    file_output.write(" ".join(f"{x:.6f}" for x in x))