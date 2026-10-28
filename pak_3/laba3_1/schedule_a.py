
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.pi / 2 - np.arctan(x)


def lagrange(x, nodes):
    result = 0
    for i in range(4):
        term = f(nodes[i])
        for j in range(4):
            if i != j:
                term *= (x - nodes[j]) / (nodes[i] - nodes[j])
        result += term
    return result


def newton(x, nodes):
    y = f(np.array(nodes, dtype=float))
    n = len(nodes)

    for j in range(1, n):
        y[j:] = (y[j:] - y[j - 1:-1]) / (
            np.array(nodes[j:]) - np.array(nodes[:-j])
        )

    result = np.zeros_like(x, dtype=float) + y[0]
    term = np.ones_like(x, dtype=float)

    for i in range(1, n):
        term *= x - nodes[i - 1]
        result += y[i] * term

    return result


nodes = [-3, -1, 1, 3]
x_star = -0.5

x = np.linspace(-3.5, 3.5, 1000)

plt.plot(x, f(x), label="Исходная функция f(x)")
plt.plot(x, lagrange(x, nodes), label="Многочлен Лагранжа")
plt.plot(x, newton(x, nodes), "--", label="Многочлен Ньютона")

plt.scatter(nodes, f(np.array(nodes)), color="red",
            label="Узлы интерполяции", zorder=3)

plt.axvline(x_star, color="gray", linestyle=":",
            label="x* = -0.5")

plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(0, color="black", linewidth=0.8)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Интерполяция функции arcctg(x)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
