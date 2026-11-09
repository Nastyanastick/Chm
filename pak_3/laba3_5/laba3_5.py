import math


with open("input.txt") as file_input:
    a, b = map(float, file_input.readline().split())
    h1, h2 = map(float, file_input.readline().split())

def f(x):
    return 1 / (x**3 + 64)


def left_rectangles(h):
    N = int(round((b - a) / h))
    result = 0

    for i in range(N):
        x = a + i * h
        result += f(x)

    return h * result


def right_rectangles(h):
    N = int(round((b - a) / h))
    result = 0

    for i in range(1, N + 1):
        x = a + i * h
        result += f(x)

    return h * result


def middle_rectangles(h):
    N = int(round((b - a) / h))
    result = 0

    for i in range(N):
        x = (a + i * h + a + (i + 1) * h) / 2
        result += f(x)

    return h * result


def trapezoids(h):
    N = int(round((b - a) / h))
    result = (f(a) + f(b)) / 2

    for i in range(1, N):
        x = a + i * h
        result += f(x)

    return h * result


def simpson(h):
    N = int(round((b - a) / h))
    result = f(a) + f(b)

    for i in range(1, N):
        x = a + i * h

        if i % 2 == 1:
            result += 4 * f(x)
        else:
            result += 2 * f(x)

    return h * result / 3


def runge_romberg(F1, F2, p):
    R = abs(F2 - F1) / ((h1 / h2)**p - 1)
    F = F2 + (F2 - F1) / ((h1 / h2)**p - 1)

    return F, R


def G(x):
    return (math.log(abs(x + 4)) / 48 - math.log(x**2 - 4 * x + 16) / 96 + math.sqrt(3) / 48 * math.atan((x - 2) / (2 * math.sqrt(3))))

res = G(b) - G(a)
# методы и их порядки точности
methods = [
    ("Левые прямоугольники", left_rectangles, 1),
    ("Правые прямоугольники", right_rectangles, 1),
    ("Средние прямоугольники", middle_rectangles, 2),
    ("Трапеции", trapezoids, 2),
    ("Симпсона", simpson, 4)
]
# вычисление результатов
results = []

for name, method, p in methods:
    F1 = method(h1)
    F2 = method(h2)

    F, R = runge_romberg(F1, F2, p)

    error1 = abs(res - F1)
    error2 = abs(res - F2)
    error_refined = abs(res - F)

    results.append((name, F1, F2, p, F, R, error1, error2, error_refined))

# поиск наиболее точного метода
best = results[0]
for item in results:
    if item[8] < best[8]:
        best = item

# вывод результатов
with open("output.txt", "w", encoding="utf-8") as file_output:
    file_output.write("1. Исходные данные\n\n")
    file_output.write(f"Границы интегрирования: [{a}; {b}]\n")
    file_output.write(f"Шаг h1 = {h1}\n")
    file_output.write(f"Шаг h2 = {h2}\n")
    file_output.write("Функция: f(x) = 1 / (x^3 + 64)\n")

    file_output.write("\n2. Точное значение интеграла\n\n")
    file_output.write(f"F = {res:.12f}\n")

    file_output.write("\n3. Результаты численного интегрирования\n\n")

    file_output.write(
        f"{'Метод':<25}"
        f"{'F(h1)':>16}"
        f"{'F(h2)':>16}"
        f"{'p':>6}"
        f"{'F уточн.':>16}"
        f"{'R':>16}\n"
    )

    for item in results:
        name, F1, F2, p, F, R = item[:6]

        file_output.write(
            f"{name:<25}"
            f"{F1:>16.10f}"
            f"{F2:>16.10f}"
            f"{p:>6}"
            f"{F:>16.10f}"
            f"{R:>16.10f}\n"
        )

    file_output.write("\n4. Абсолютные погрешности\n\n")
    file_output.write(
        f"{'Метод':<25}"
        f"{'h1':>16}"
        f"{'h2':>16}"
        f"{'Рунге-Ромберг':>20}\n"
    )

    for item in results:
        name, error1, error2, error_refined = (
            item[0], item[6], item[7], item[8]
        )

        file_output.write(
            f"{name:<25}"
            f"{error1:>16.10f}"
            f"{error2:>16.10f}"
            f"{error_refined:>20.10f}\n"
        )

    file_output.write("\n5. Наиболее точный метод\n\n")
    file_output.write(f"Метод: {best[0]}\n")
    file_output.write(f"Уточнённое значение: {best[4]:.12f}\n")
    file_output.write(f"Абсолютная погрешность: {best[8]:.12f}\n")
