with open("input.txt") as file_input:
    x = list(map(float, file_input.readline().split()))
    y = list(map(float, file_input.readline().split()))
    X = float(file_input.readline())

N = len(x) - 1

# поиск узла, совпадающего с X
index = -1

for i in range(N + 1):
    if abs(x[i] - X) < 1e-9:
        index = i
        break

def derivatives(i):
    h = x[i + 1] - x[i]

    first = (y[i + 1] - y[i]) / h

    second = (
        y[i + 2] - 2 * y[i + 1] + y[i]
    ) / (h**2)

    first_derivative = first + (
        (y[i + 2] - 2 * y[i + 1] + y[i]) / (2 * h**2)
    ) * (2 * X - x[i] - x[i + 1])

    return first_derivative, second


left_derivative = (y[index] - y[index - 1]) / (x[index] - x[index - 1])

right_derivative = (y[index + 1] - y[index]) / (x[index + 1] - x[index])

first_derivative, second_derivative = derivatives(index - 1)

with open("output.txt", "w", encoding="utf-8") as file_output:
    file_output.write("1. Исходная таблица значений\n\n")

    for i in range(N + 1):
        file_output.write(
            f"x[{i}] = {x[i]:.4f}, y[{i}] = {y[i]:.4f}\n"
        )

    file_output.write(f"\n2. Точка вычисления: X* = {X:.4f}\n")

    file_output.write("\n3. Узлы, использованные при вычислении\n\n")
    file_output.write(
        f"Для левой производной: x[{index - 1}], x[{index}]\n"
    )
    file_output.write(
        f"Для правой производной: x[{index}], x[{index + 1}]\n"
    )
    file_output.write(
        f"Для трёхточечных формул: "
        f"x[{index - 1}], x[{index}], x[{index + 1}]\n"
    )

    file_output.write("\n4. Первая производная\n\n")
    file_output.write(
        f"Левая: {left_derivative:.8f}\n"
    )
    file_output.write(
        f"Правая: {right_derivative:.8f}\n"
    )
    file_output.write(
        f"По трём узлам: {first_derivative:.8f}\n"
    )

    file_output.write("\n5. Вторая производная\n\n")
    file_output.write(
        f"По трём узлам: {second_derivative:.8f}\n"
    )