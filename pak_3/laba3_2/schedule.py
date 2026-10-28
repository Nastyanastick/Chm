import numpy as np
import matplotlib.pyplot as plt

nodes = [-3.0, -1.0, 1.0, 3.0, 5.0]
values = [2.8198, 2.3562, 0.78540, 0.32175, 0.19740]
x_star = -0.5

a = [2.8198000000, 2.3562000000, 0.7854000000, 0.3217500000]
b = [-0.0470026786, -0.6013946429, -0.5990187500, -0.0542053571]
c = [0.0, -0.2771959821, 0.2783839286, -0.0059772321, 0.0]
d = [-0.0461993304, 0.0925966518, -0.0473935268, 0.0009962054]

plt.figure(figsize=(10, 6))

for i in range(len(nodes) - 1):
    x = np.linspace(nodes[i], nodes[i + 1], 200)

    y = (
        a[i]
        + b[i] * (x - nodes[i])
        + c[i] * (x - nodes[i])**2
        + d[i] * (x - nodes[i])**3
    )

    plt.plot(x, y, label=f"Сплайн S{i}(x)")

plt.scatter(
    nodes,
    values,
    color="red",
    label="Исходные узлы",
    zorder=3
)

interval = 1
s_star = (
    a[interval]
    + b[interval] * (x_star - nodes[interval])
    + c[interval] * (x_star - nodes[interval])**2
    + d[interval] * (x_star - nodes[interval])**3
)

plt.scatter(
    x_star,
    s_star,
    color="green",
    label=f"(x*, S(x*)) = ({x_star}, {s_star:.6f})",
    zorder=4
)

plt.axvline(x_star, color="gray", linestyle=":")

plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(0, color="black", linewidth=0.8)

plt.xlabel("x")
plt.ylabel("S(x)")
plt.title("Кубический сплайн")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()