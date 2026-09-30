import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return np.exp(x) - x**3 + 3 * x**2 - 2 * x - 3


def y1(x):
    return np.exp(x)


def y2(x):
    return x**3 - 3 * x**2 + 2 * x + 3


x = np.linspace(0, 1.5, 10000)      # массив из 500 чисел
plt.plot(x, y1(x), label="y = e^x")
plt.plot(x, y2(x), label="y = x**3 - 3 * x**2 + 2 * x + 3")

plt.axhline(0, color="black", linewidth=0.8)    # горизонтальная ось
plt.axvline(0, color="black", linewidth=0.8)    # вертикальная ось

plt.xlabel("x")     # подписывает x
plt.ylabel("y")
plt.title("Графическое решение уравнения")
plt.grid(True)      # включает сетку
plt.legend()
plt.show()      # отображает график
