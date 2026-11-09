
import math


a = 1.0
b = 2.0
h = 0.1
delta = 1e-8

def f(x, y, z):
    return z


def g(x, y, z):
    return (x * z - y) / (x * (x - 1))


def get_initial_values(s, eps, terms=12):   # 12 членов ряда начально значение выбрали такое
    coefficients = [0.0] * (terms + 1)

    coefficients[0] = 2.0
    coefficients[1] = 2.0
    coefficients[2] = s / 2

    for n in range(2, terms):
        coefficients[n + 1] = -(n - 1) / (n + 1) * coefficients[n]

    y = 0.0
    z = 0.0

    for n in range(terms + 1):
        y += coefficients[n] * eps**n

        if n > 0:
            z += n * coefficients[n] * eps**(n - 1)

    return y, z


def runge_kutta(h, s, eps, terms=12):
    x0 = a + eps
    N = round((b - x0) / h)
    step = (b - x0) / N

    x = [x0 + k * step for k in range(N + 1)]

    y = [0.0] * (N + 1)
    z = [0.0] * (N + 1)

    y[0], z[0] = get_initial_values(s, eps, terms)

    for k in range(N):
        K1 = step * f(x[k], y[k], z[k])
        L1 = step * g(x[k], y[k], z[k])

        K2 = step * f(x[k] + step / 2, y[k] + K1 / 2, z[k] + L1 / 2)
        L2 = step * g(x[k] + step / 2, y[k] + K1 / 2, z[k] + L1 / 2)

        K3 = step * f(x[k] + step / 2, y[k] + K2 / 2, z[k] + L2 / 2)
        L3 = step * g(x[k] + step / 2, y[k] + K2 / 2, z[k] + L2 / 2)

        K4 = step * f(x[k] + step, y[k] + K3, z[k] + L3)
        L4 = step * g(x[k] + step, y[k] + K3, z[k] + L3)

        y[k + 1] = y[k] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        z[k + 1] = z[k] + (L1 + 2 * L2 + 2 * L3 + L4) / 6

    return x, y, z


def calculate_residual(s, h, eps, terms=12):    # наша Ф
    x, y, z = runge_kutta(h, s, eps, terms)

    return 2 * z[-1] - y[-1] - 1


def shooting_method(h, eps, delta, terms=12):
    s0 = 0.0
    s1 = 2.0

    phi0 = calculate_residual(s0, h, eps, terms)
    phi1 = calculate_residual(s1, h, eps, terms)

    iterations = 0
    max_iterations = 100

    while abs(phi1) > delta and iterations < max_iterations:
        denominator = phi1 - phi0

        if abs(denominator) < 1e-14:
            raise ValueError("Невозможно продолжить метод секущих")

        s2 = s1 - phi1 * (s1 - s0) / denominator

        s0 = s1
        phi0 = phi1

        s1 = s2
        phi1 = calculate_residual(s1, h, eps, terms)

        iterations += 1

    if abs(phi1) > delta:
        raise ValueError("Метод секущих не сошелся")

    x, y, z = runge_kutta(h, s1, eps, terms)

    return s1, x, y, z, phi1, iterations


def gauss(A, D):
    n = len(D)
    A = [row[:] for row in A]
    D = D[:]

    for k in range(n):
        max_row = k

        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[max_row][k]):
                max_row = i

        A[k], A[max_row] = A[max_row], A[k]
        D[k], D[max_row] = D[max_row], D[k]

        if abs(A[k][k]) < 1e-14:
            raise ValueError("Система не имеет единственного решения")

        for i in range(k + 1, n):
            factor = A[i][k] / A[k][k]

            for j in range(k, n):
                A[i][j] -= factor * A[k][j]

            D[i] -= factor * D[k]

    y = [0.0] * n

    for i in range(n - 1, -1, -1):
        y[i] = D[i]

        for j in range(i + 1, n):
            y[i] -= A[i][j] * y[j]

        y[i] /= A[i][i]

    return y


def finite_difference(h):
    N = round((b - a) / h)
    h = (b - a) / N

    x = [a + i * h for i in range(N + 1)]

    A = [[0.0] * N for _ in range(N)]
    D = [0.0] * N

    y0 = 2.0

    for i in range(1, N):
        xi = x[i]

        c1 = xi * (xi - 1) / h**2 + xi / (2 * h)
        c2 = -2 * xi * (xi - 1) / h**2 + 1
        c3 = xi * (xi - 1) / h**2 - xi / (2 * h)

        row = i - 1

        if i == 1:
            D[row] = -c1 * y0
        else:
            A[row][row - 1] = c1

        A[row][row] = c2
        A[row][row + 1] = c3

    A[N - 1][N - 3] = 1 / h
    A[N - 1][N - 2] = -4 / h
    A[N - 1][N - 1] = 3 / h - 1
    D[N - 1] = 1.0

    solution = gauss(A, D)

    y = [y0] + solution

    return x, y


def exact(x):
    return 1 + x + x * math.log(x)


eps = 0.01
terms = 12

s_h, x_h, y_h, z_h, residual_h, iterations_h = shooting_method(h, eps, delta, terms)
s_half, x_half, y_half, z_half, residual_half, iterations_half = shooting_method(h / 2, eps, delta, terms)

x_fd_h, y_fd_h = finite_difference(h)
x_fd_half, y_fd_half = finite_difference(h / 2)

residual_fd_h = 2 * (3 * y_fd_h[-1] - 4 * y_fd_h[-2] + y_fd_h[-3]) / (2 * h) - y_fd_h[-1] - 1
residual_fd_half = 2 * (3 * y_fd_half[-1] - 4 * y_fd_half[-2] + y_fd_half[-3]) / h - y_fd_half[-1] - 1


with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("МЕТОД СТРЕЛЬБЫ\n")

    file_output.write(f"Шаг h = {h}\n")
    file_output.write(f"Найденный параметр s = {s_h:.12f}\n\n")
    file_output.write(f"Количество итераций = {iterations_h}\n")
    file_output.write(f"Невязка = {residual_h:.12f}\n\n")

    file_output.write(f"Шаг h/2 = {h / 2}\n")
    file_output.write(f"Найденный параметр s = {s_half:.12f}\n")
    file_output.write(f"Количество итераций = {iterations_half}\n")
    file_output.write(f"Невязка = {residual_half:.12f}\n\n")

    file_output.write(
        f"{'x':>10} {'y(h)':>16} {'y(h/2)':>16} "
        f"{'y точное':>16} {'Погр. Р-Р':>16} {'Погр. точная':>16}\n"
    )

    for k in range(len(x_h)):
        j = 2 * k

        y_exact = exact(x_h[k])
        rr_error = abs((y_half[j] - y_h[k]) / 15)
        exact_error = abs(y_half[j] - y_exact)

        file_output.write(
            f"{x_h[k]:10.4f} "
            f"{y_h[k]:16.10f} "
            f"{y_half[j]:16.10f} "
            f"{y_exact:16.10f} "
            f"{rr_error:16.10f} "
            f"{exact_error:16.10f}\n"
        )

    file_output.write("\n\nМЕТОД КОНЕЧНЫХ РАЗНОСТЕЙ\n\n")

    file_output.write(f"Шаг h = {h}\n")
    file_output.write(f"Невязка = {residual_fd_h:.12f}\n\n")

    file_output.write(f"Шаг h/2 = {h / 2}\n")
    file_output.write(f"Невязка = {residual_fd_half:.12f}\n\n")

    file_output.write(
        f"{'x':>10} {'y(h)':>16} {'y(h/2)':>16} "
        f"{'y точное':>16} {'Погр. Р-Р':>16} {'Погр. точная':>16}\n"
    )

    for k in range(len(x_fd_h)):
        j = 2 * k

        y_exact = exact(x_fd_h[k])
        rr_error = abs((y_fd_half[j] - y_fd_h[k]) / 3)
        exact_error = abs(y_fd_half[j] - y_exact)

        file_output.write(
            f"{x_fd_h[k]:10.4f} "
            f"{y_fd_h[k]:16.10f} "
            f"{y_fd_half[j]:16.10f} "
            f"{y_exact:16.10f} "
            f"{rr_error:16.10f} "
            f"{exact_error:16.10f}\n"
        )
