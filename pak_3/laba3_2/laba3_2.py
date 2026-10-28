import math

with open ("input.txt") as file_input:
    x = list(map(float, file_input.readline().split()))
    y = list(map(float, file_input.readline().split()))
    x_star = float(file_input.readline())

n = len(x) - 1
h = [0.0] * n

for i in range(n):
    h[i] = x[i + 1] - x[i]

delta = [0.0] * n

for i in range(n):
    delta[i] = (y[i + 1] - y[i]) / h[i]

A = [[0.0] * 3 for i in range(3)]
B = [0.0] * 3

# внутр узлы
for i in range(1, n):
    A[i - 1][i - 1] = 2 * (h[i - 1] + h[i])

    if i > 1:
        A[i - 1][i - 2] = h[i - 1]

    if i < n - 1:
        A[i - 1][i] = h[i]

    B[i - 1] = 3 * (delta[i] - delta[i - 1])

# метод гаусса для решения системы
for i in range(3):
    for j in range(i + 1, 3):
        factor = A[j][i] / A[i][i]

        for k in range(i, 3):
            A[j][k] -= factor * A[i][k]

        B[j] -= factor * B[i]

c = [0.0] * 5

for i in range(2, -1, -1):
    c[i + 1] = B[i]

    for j in range(i + 1, 3):
        c[i + 1] -= A[i][j] * c[j + 1]

    c[i + 1] /= A[i][i]

a = [0.0] * n
b = [0.0] * n
d = [0.0] * n

for i in range(n):
    a[i] = y[i]
    d[i] = (c[i + 1] - c[i]) / (3 * h[i])
    b[i] = delta[i] - c[i] * h[i] - d[i] * h[i]**2

# промежуток где x*
interval = 0
for i in range(n):
    if x[i] <= x_star <= x[i + 1]:
        interval = i + 1
        s = (
            a[i]
            + b[i] * (x_star - x[i])
            + c[i] * (x_star - x[i])**2
            + d[i] * (x_star - x[i])**3
        )
        break

with open("output.txt", "w", encoding="utf-8") as file_output:
    file_output.write("1. Исходные узлы и значения функции\n\n")

    for i in range(len(x)):
        file_output.write(f"x[{i}] = {x[i]}, y[{i}] = {y[i]}\n")

    file_output.write("\n2-3. Коэффициенты и интервалы сплайна\n\n")

    for i in range(n):
        file_output.write(f"Участок {i + 1}: [{x[i]}, {x[i + 1]}]\n")
        file_output.write(f"a[{i}] = {a[i]:.10f}\n")
        file_output.write(f"b[{i}] = {b[i]:.10f}\n")
        file_output.write(f"c[{i}] = {c[i]:.10f}\n")
        file_output.write(f"d[{i}] = {d[i]:.10f}\n")
        file_output.write(
            f"S{i}(x) = {a[i]:.10f} "
            f"+ ({b[i]:.10f}) * (x - ({x[i]})) "
            f"+ ({c[i]:.10f}) * (x - ({x[i]}))^2 "
            f"+ ({d[i]:.10f}) * (x - ({x[i]}))^3\n\n"
        )

    file_output.write("4. Значение сплайна в точке x*\n")
    file_output.write(f"x* = {x_star}\n")
    file_output.write(f"S(x*) = {s:.10f}\n")

    file_output.write("\n5. Номер интервала\n")
    file_output.write(
        f"Интервал номер {interval}: "
        f"[{x[interval - 1]}, {x[interval]}]\n"
    )

    file_output.write("\n6. Проверка прохождения сплайна через узлы\n\n")

    for i in range(n):
        left_value = a[i]
        right_value = (
            a[i] + b[i] * h[i]
            + c[i] * h[i]**2 + d[i] * h[i]**3
        )

        file_output.write(
            f"Узел x[{i}]: S(x[{i}]) = {left_value:.10f}, "
            f"y[{i}] = {y[i]:.10f}\n"
        )

        if i == n - 1:
            file_output.write(
                f"Узел x[{i + 1}]: S(x[{i + 1}]) = "
                f"{right_value:.10f}, y[{i + 1}] = "
                f"{y[i + 1]:.10f}\n"
            )

    file_output.write("\n7. Проверка граничных условий\n")

    second_left = 2 * c[0]
    second_right = 2 * c[n - 1] + 6 * d[n - 1] * h[n - 1]

    file_output.write(f"S''(x[0]) = {second_left:.10f}\n")
    file_output.write(f"S''(x[{n}]) = {second_right:.10f}\n")