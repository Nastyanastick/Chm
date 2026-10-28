import math

with open("input.txt") as file_input:
    xa = list(map(int, file_input.readline().split()))
    xb = list(map(int, file_input.readline().split()))
    x = float(file_input.readline())
    fxa = [0.0] * 4
    fxb = [0.0] * 4

La = [0.0] * 4
Lb = [0.0] * 4

for i in range(4):
    fxa[i] = math.pi / 2 - math.atan(xa[i])
    fxb[i] = math.pi / 2 - math.atan(xb[i])

# лагранж
for i in range(4):
    la = [1.0]
    lb = [1.0]
    for j in range(4):
        if i != j:
            # умножаем текущий многочлен на (x - xj)
            new_la = [0.0] * (len(la) + 1)
            new_lb = [0.0] * (len(lb) + 1)

            for k in range(len(la)):
                new_la[k] -= xa[j] * la[k]
                new_la[k + 1] += la[k]

            for k in range(len(lb)):
                new_lb[k] -= xb[j] * lb[k]
                new_lb[k + 1] += lb[k]

            la = new_la

            lb = new_lb
    # делим на произведение (xi - xj)
    den_a = 1.0
    den_b = 1.0
    for j in range(4):
        if i != j:
            den_a *= xa[i] - xa[j]
            den_b *= xb[i] - xb[j]
    for k in range(len(la)):
        la[k] /= den_a
    for k in range(len(lb)):
        lb[k] /= den_b

    for k in range(4):
        La[k] += fxa[i] * la[k]
        Lb[k] += fxb[i] * lb[k]

# таблицы разделяемых разностей
dda = [[0.0] * 4 for i in range(4)]
ddb = [[0.0] * 4 for i in range(4)]
for i in range(4):
    dda[i][0] = fxa[i]
    ddb[i][0] = fxb[i]
# вычисляем разности
for j in range(1, 4):   # порядок разности
    for i in range(4 - j):  # номер строки
        dda[i][j] = (dda[i][j - 1] - dda[i + 1][j - 1]) / (xa[i] - xa[i + j])
        ddb[i][j] = (ddb[i][j - 1] - ddb[i + 1][j - 1]) / (xb[i] - xb[i + j])
Na = dda[0][:]
Nb = ddb[0][:]

# вывод для нашего x
Na_x = 0.0
Nb_x = 0.0

for i in range(4):
    pa = 1.0
    pb = 1.0

    for j in range(i):
        pa *= x - xa[j]
        pb *= x - xb[j]

    Na_x += dda[0][i] * pa
    Nb_x += ddb[0][i] * pb


La_x = sum(La[k] * x**k for k in range(4))
Lb_x = sum(Lb[k] * x**k for k in range(4))

# погрешность
def fourth_derivative(x):
    return 24 * x * (x**2 - 1) / (1 + x**2)**4

critical_points = [
    -3.0,
    3.0,
    -math.sqrt((5 - 2 * math.sqrt(5)) / 5),
    math.sqrt((5 - 2 * math.sqrt(5)) / 5),
    -math.sqrt((5 + 2 * math.sqrt(5)) / 5),
    math.sqrt((5 + 2 * math.sqrt(5)) / 5)
]

M4 = max(abs(fourth_derivative(t)) for t in critical_points)

omega_a = 1.0
omega_b = 1.0

for i in range(4):
    omega_a *= x - xa[i]
    omega_b *= x - xb[i]

error_bound_a = M4 / math.factorial(4) * abs(omega_a)
error_bound_b = M4 / math.factorial(4) * abs(omega_b)

fx = math.pi / 2 - math.atan(x)
actual_error_a = abs(fx - La_x)
actual_error_b = abs(fx - Lb_x)

with open("output.txt", "w", encoding="utf-8") as file_output:

    file_output.write("1. Исходные данные\n\n")
    file_output.write("Вариант а:\n")
    for i in range(4):
        file_output.write(f"x[{i}] = {xa[i]}, y[{i}] = {fxa[i]:.6f}\n")
    file_output.write("\nВариант б:\n")
    for i in range(4):
        file_output.write(f"x[{i}] = {xb[i]}, y[{i}] = {fxb[i]:.6f}\n")

    # 2. многочлен Лагранжа
    file_output.write("\n2. Многочлены Лагранжа\n\n")
    file_output.write("Вариант а:\n")
    file_output.write(
        f"L(x) = {La[0]:.6f} "
        f"+ ({La[1]:.6f}) * x "
        f"+ ({La[2]:.6f}) * x^2 "
        f"+ ({La[3]:.6f}) * x^3\n"
    )

    file_output.write("\nВариант б:\n")
    file_output.write(
        f"L(x) = {Lb[0]:.6f} "
        f"+ ({Lb[1]:.6f}) * x "
        f"+ ({Lb[2]:.6f}) * x^2 "
        f"+ ({Lb[3]:.6f}) * x^3\n"
    )

    # 3. Таблицы разностей
    file_output.write("\n3. Таблица разделённых разностей для варианта а:\n")
    for i in range(4):
        file_output.write(" ".join(f"{dda[i][j]:.6f}" for j in range(4 - i)) + "\n")

    file_output.write("\nТаблица разделённых разностей для варианта б:\n")
    for i in range(4):
        file_output.write(" ".join(f"{ddb[i][j]:.6f}" for j in range(4 - i)) + "\n")

    # 4. Многочлены Ньютона
    file_output.write("\n4. Многочлены Ньютона\n\n")

    file_output.write("Вариант а:\n")
    file_output.write(
        f"N(x) = {Na[0]:.6f} "
        f"+ ({Na[1]:.6f}) * (x - ({xa[0]})) "
        f"+ ({Na[2]:.6f}) * (x - ({xa[0]})) * (x - ({xa[1]})) "
        f"+ ({Na[3]:.6f}) * (x - ({xa[0]})) * (x - ({xa[1]})) * (x - ({xa[2]}))\n"
    )

    file_output.write("\nВариант б:\n")
    file_output.write(
        f"N(x) = {Nb[0]:.6f} "
        f"+ ({Nb[1]:.6f}) * (x - ({xb[0]})) "
        f"+ ({Nb[2]:.6f}) * (x - ({xb[0]})) * (x - ({xb[1]})) "
        f"+ ({Nb[3]:.6f}) * (x - ({xb[0]})) * (x - ({xb[1]})) * (x - ({xb[2]}))\n"
    )

    file_output.write("\n5. Значение многочлена Лагранжа в точке x*:\n")
    file_output.write(f"Вариант а: L(x*) = {La_x:.10f}\n")
    file_output.write(f"Вариант б: L(x*) = {Lb_x:.10f}\n")

    file_output.write("\n6. Значение многочлена Ньютона в точке x*:\n")
    file_output.write(f"Вариант а: N(x*) = {Na_x:.10f}\n")
    file_output.write(f"Вариант б: N(x*) = {Nb_x:.10f}\n")

    file_output.write("\n7. Точное значение функции:\n")
    fx = math.pi / 2 - math.atan(x)
    file_output.write(f"f(x*) = {fx:.10f}\n")

    file_output.write("\n8. Оценка погрешности:\n")
    file_output.write(f"M4 = {M4:.6f}\n")

    file_output.write("\nВариант а:\n")
    file_output.write(f"|omega4(x*)| = {abs(omega_a):.6f}\n")
    file_output.write(f"Оценка по формуле = {error_bound_a:.10f}\n")
    file_output.write(
        f"|f(x*) - L(x*)| = {actual_error_a:.10f} "
        f"<= {error_bound_a:.10f}\n"
    )

    file_output.write("\nВариант б:\n")
    file_output.write(f"|omega4(x*)| = {abs(omega_b):.6f}\n")
    file_output.write(f"Оценка по формуле = {error_bound_b:.10f}\n")
    file_output.write(
        f"|f(x*) - L(x*)| = {actual_error_b:.10f} "
        f"<= {error_bound_b:.10f}\n"
    )

    file_output.write("\n9. Разность результатов методов:\n")
    file_output.write(f"Вариант а: |L(x*) - N(x*)| = {abs(La_x - Na_x):.10f}\n")
    file_output.write(f"Вариант б: |L(x*) - N(x*)| = {abs(Lb_x - Nb_x):.10f}\n")