import math


h = 0.1
a = 0.0
b = 1.0
y0 = 2.0
z0 = 0.0

def f(x, y, z):
    return z


def g(x, y, z):
    return y * math.cos(x)**2 - z * math.tan(x)


def exact(x):
    return math.exp(math.sin(x)) + math.exp(-math.sin(x))

def exact_derivative(x):
    return math.cos(x) * (math.exp(math.sin(x)) - math.exp(-math.sin(x)))


def euler(h):
    N = round((b - a) / h)
    x = [a + k * h for k in range(N + 1)]

    y = [0.0] * (N + 1)
    z = [0.0] * (N + 1)

    y[0] = y0
    z[0] = z0

    for k in range(N):
        y[k + 1] = y[k] + h * f(x[k], y[k], z[k])
        z[k + 1] = z[k] + h * g(x[k], y[k], z[k])

    return x, y, z



def runge_kutta(h):
    N = round((b - a) / h)

    x = [a + k * h for k in range(N + 1)]

    y = [0.0] * (N + 1)
    z = [0.0] * (N + 1)

    y[0] = y0
    z[0] = z0

    for k in range(N):
        K1 = h * f(x[k], y[k], z[k])
        L1 = h * g(x[k], y[k], z[k])

        K2 = h * f(x[k] + h / 2, y[k] + K1 / 2, z[k] + L1 / 2)
        L2 = h * g(x[k] + h / 2, y[k] + K1 / 2, z[k] + L1 / 2)

        K3 = h * f(x[k] + h / 2, y[k] + K2 / 2, z[k] + L2 / 2)
        L3 = h * g(x[k] + h / 2, y[k] + K2 / 2, z[k] + L2 / 2)

        K4 = h * f(x[k] + h, y[k] + K3, z[k] + L3)
        L4 = h * g(x[k] + h, y[k] + K3, z[k] + L3)

        y[k + 1] = y[k] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        z[k + 1] = z[k] + (L1 + 2 * L2 + 2 * L3 + L4) / 6

    return x, y, z


def adams(h):
    N = round((b - a) / h)

    x, y, z = runge_kutta(h)

    for k in range(3, N):
        y[k + 1] = y[k] + h / 24 * (
            55 * f(x[k], y[k], z[k])
            - 59 * f(x[k - 1], y[k - 1], z[k - 1])
            + 37 * f(x[k - 2], y[k - 2], z[k - 2])
            - 9 * f(x[k - 3], y[k - 3], z[k - 3])
        )

        z[k + 1] = z[k] + h / 24 * (
            55 * g(x[k], y[k], z[k])
            - 59 * g(x[k - 1], y[k - 1], z[k - 1])
            + 37 * g(x[k - 2], y[k - 2], z[k - 2])
            - 9 * g(x[k - 3], y[k - 3], z[k - 3])
        )

    return x, y, z


def runge_romberg(y_h, z_h, y_h2, z_h2, p):
    error_y = []
    error_z = []

    for k in range(len(y_h)):
        error_y.append(abs(y_h2[2 * k] - y_h[k]) / (2**p - 1))
        error_z.append(abs(z_h2[2 * k] - z_h[k]) / (2**p - 1))

    return error_y, error_z

x_e, y_e, z_e = euler(h)
x_e2, y_e2, z_e2 = euler(h / 2)

x_rk, y_rk, z_rk = runge_kutta(h)
x_rk2, y_rk2, z_rk2 = runge_kutta(h / 2)

x_ad, y_ad, z_ad = adams(h)
x_ad2, y_ad2, z_ad2 = adams(h / 2)


error_e_y, error_e_z = runge_romberg(y_e, z_e, y_e2, z_e2, 1)
error_rk_y, error_rk_z = runge_romberg(y_rk, z_rk, y_rk2, z_rk2, 4)
error_ad_y, error_ad_z = runge_romberg(y_ad, z_ad, y_ad2, z_ad2, 4)

with open("output.txt", "w", encoding="utf-8") as file_output:

    methods = [
        ("Метод Эйлера", x_e, y_e, z_e, error_e_y, error_e_z),
        ("Метод Рунге — Кутты 4-го порядка", x_rk, y_rk, z_rk, error_rk_y, error_rk_z),
        ("Метод Адамса 4-го порядка", x_ad, y_ad, z_ad, error_ad_y, error_ad_z)
    ]

    for name, x, y, z, error_y, error_z in methods:
        file_output.write(f"\n{name}\n")
        file_output.write("-" * 110 + "\n")

        file_output.write(
            f"{'x':>8} {'y прибл.':>14} {'y точное':>14} "
            f"{'Ошибка y':>14} {'Рунге y':>14} "
            f"{'z прибл.':>14} {'z точное':>14} "
            f"{'Ошибка z':>14} {'Рунге z':>14}\n"
        )

        for k in range(len(x)):
            y_exact = exact(x[k])
            z_exact = exact_derivative(x[k])

            file_output.write(
                f"{x[k]:8.2f} "
                f"{y[k]:14.8f} "
                f"{y_exact:14.8f} "
                f"{abs(y_exact - y[k]):14.12f} "
                f"{error_y[k]:14.12f} "
                f"{z[k]:14.8f} "
                f"{z_exact:14.8f} "
                f"{abs(z_exact - z[k]):14.12f} "
                f"{error_z[k]:14.12f}\n"
            )