import numpy as np
import matplotlib.pyplot as plt


def y1(x2):
    return 2 + np.cos(x2)


def y2(x1):
    return 2 + np.sin(x1)


x1 = np.linspace(0, 3.5, 10000)
x2 = np.linspace(0, 3.5, 10000)

plt.plot(y1(x2), x2, label="x1 = 2 + cos(x2)")
plt.plot(x1, y2(x1), label="x2 = 2 + sin(x1)")

plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(0, color="black", linewidth=0.8)

plt.xlabel("x1")
plt.ylabel("x2")
plt.title("Графическое решение системы нелинейных уравнений")

plt.grid(True)
plt.legend()
plt.show()
