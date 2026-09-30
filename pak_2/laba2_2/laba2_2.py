import math


def f1(x1, x2):
    return x1 - math.cos(x2) - 2


def f2(x1, x2):
    return x2 - math.sin(x1) - 2


def det(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


with open("input.txt") as file_input:
    x1_n = float(file_input.readline())
    x2_n = float(file_input.readline())
    x1_iter = x1_n
    x2_iter = x2_n
    eps = float(file_input.readline())


# ньютон
k_n = 0

while True:

    F1 = f1(x1_n, x2_n)
    F2 = f2(x1_n, x2_n)

    J = [
        [1, math.sin(x2_n)],
        [-math.cos(x1_n), 1]
    ]

    A1 = [
        [F1, math.sin(x2_n)],
        [F2, 1]
    ]

    A2 = [
        [1, F1],
        [-math.cos(x1_n), F2]
    ]

    det_J = det(J)

    if abs(det_J) < 1e-12:
        raise ValueError("Определитель матрицы Якоби равен нулю. Метод Ньютона не может продолжить вычисления.")

    new_x1_n = x1_n - det(A1) / det_J
    new_x2_n = x2_n - det(A2) / det_J

    delta_n = max(abs(new_x1_n - x1_n), abs(new_x2_n - x2_n))

    k_n += 1

    x1_n = new_x1_n
    x2_n = new_x2_n

    if delta_n <= eps:
        break


# метод простой итерации
k_iter = 0
q = math.cos(1.05)

while True:

    new_x1_iter = 2 + math.cos(x2_iter)
    new_x2_iter = 2 + math.sin(x1_iter)

    delta_iter = max(
        abs(new_x1_iter - x1_iter),
        abs(new_x2_iter - x2_iter)
    )

    k_iter += 1

    x1_iter = new_x1_iter
    x2_iter = new_x2_iter

    if q / (1 - q) * delta_iter <= eps:
        break


with open("output.txt", "w", encoding="utf-8") as file_output:
    file_output.write(f"Корень x1 (метод Ньютона): {x1_n:.6f}\n")
    file_output.write(f"Корень x2 (метод Ньютона): {x2_n:.6f}\n")
    file_output.write(f"Количество итераций (метод Ньютона): {k_n}\n")

    file_output.write(f"Корень x1 (метод простой итерации): {x1_iter:.6f}\n")
    file_output.write(f"Корень x2 (метод простой итерации): {x2_iter:.6f}\n")
    file_output.write(f"Количество итераций (метод простой итерации): {k_iter}\n")