import math


def f(x):
    return math.exp(x) - x**3 + 3 * x**2 - 2 * x - 3


def phi(x):
    return x - 0.25 * f(x)


def df(x):
    return math.exp(x) - 3 * x**2 + 6 * x - 2


with open("input.txt") as file_input:
    x0_iter = float(file_input.readline())
    x0_n = float(file_input.readline())
    eps = float(file_input.readline())


iteration_iter = 0
iteration_n = 0
x = x0_iter
q = 0.0705

while True:
    new_x_iter = phi(x)
    iteration_iter += 1

    if q / (1 - q) * abs(new_x_iter - x) < eps:
        break

    x = new_x_iter


x = x0_n
while True:
    new_x_n = x - f(x) / df(x)
    iteration_n += 1

    if abs(new_x_n - x) < eps:
        break

    x = new_x_n


with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write(f"Корень (метод простых итераций): {new_x_iter:.6f}\n")
    file_output.write(f"Количество итераций (метод простых итераций): {iteration_iter}\n")
    file_output.write(f"Корень (метод Ньютона): {new_x_n:.6f}\n")
    file_output.write(f"Количество итераций (метод Ньютона): {iteration_n}\n")