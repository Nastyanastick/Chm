import math
import matplotlib.pyplot as plt

import laba4_2 as calc


x_exact = [1 + i * 0.01 for i in range(101)]
y_exact = [1 + x + x * math.log(x) for x in x_exact]


plt.figure(figsize=(10, 6))

plt.plot(x_exact, y_exact, label="Точное решение")

plt.plot(calc.x_h, calc.y_h, "o", markersize=7, label="Метод стрельбы, шаг 0.1")

plt.plot(calc.x_half, calc.y_half, "o", markersize=5, label="Метод стрельбы, шаг 0.05")

plt.plot(calc.x_fd_h, calc.y_fd_h, "o", markersize=7, label="Конечно-разностный метод, шаг 0.1")

plt.plot(calc.x_fd_half, calc.y_fd_half, "o", markersize=5, label="Конечно-разностный метод, шаг 0.05")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Сравнение численных методов с точным решением")

plt.grid()
plt.legend()
plt.tight_layout()
plt.show()