import math
import matplotlib.pyplot as plt


h = 0.05
a = 0.0
b = 1.0
y0 = 2.0
z0 = 0.0

N = round((b - a) / h)


def f(x, y, z):
    return z


def g(x, y, z):
    return y * math.cos(x)**2 - z * math.tan(x)


def exact(x):
    return math.exp(math.sin(x)) + math.exp(-math.sin(x))


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


x_e, y_e, z_e = euler(h)
x_rk, y_rk, z_rk = runge_kutta(h)
x_ad, y_ad, z_ad = adams(h)


x_exact = [a + k * 0.001 for k in range(1001)]
y_exact = [exact(x) for x in x_exact]


plt.figure(figsize=(10, 6))

plt.plot(x_e, y_e, 'o', color='blue', label="Метод Эйлера")
plt.plot(x_e, y_e, color='blue', linewidth=1.5, label="Эйлер, h = 0.05")

plt.plot(x_rk, y_rk, 'o', color='orange',
         label="Метод Рунге-Кутты 4-го порядка")
plt.plot(x_rk, y_rk, color='orange', linewidth=1.5,
         label="Рунге-Кутта, h = 0.05")

plt.plot(x_ad, y_ad, 'o', color='green',
         label="Метод Адамса 4-го порядка")
plt.plot(x_ad, y_ad, color='green', linestyle='--',
         linewidth=1.5, label="Адамс, h = 0.05")

plt.plot(x_exact, y_exact, color='red',
         linewidth=1.8, label="Точное решение")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Сравнение численных методов с точным решением")

plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()
