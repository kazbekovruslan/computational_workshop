import numpy as np
from tabulate import tabulate
from math import factorial, tan, sqrt
import math


# Точное решение (вариант 13)
def exact_solution(x):
    sqrt_7 = sqrt(7)
    return 0.25 + sqrt_7 / 4 * tan(sqrt_7 / 2 * x)


# def exact_solution(x):
#     return (math.sin(x) - math.cos(x) + 1 / math.exp(x)) / 2
# x0, y0 = 0, 0

# Функция, определяющая правую часть уравнения
def f(x, y):
    return -y + 2 * y**2 + 1
# x0, y0 = 0, 0.25


# def f(x, y):
#     return math.sin(x) - y
# # x0, y0 = 0, 0


def derivative(i):
    arr = [0.25, 0.875, -0.4375, 0.984375, -2.7890625, -3.7734375]
    return arr[i]


# def derivative(i):
#     arr = [0, 0, 1, -1, 0, 0]
#     return arr[i]


# Метод Тейлора
def taylor_method(x0, y0, h, N):
    x_vals = [x0 + k * h for k in range(-2, N + 1)]
    y_vals = [0] * len(x_vals)
    y_vals[2] = y0
    for k in range(-2, N+1):
        x_k = x0 + k * h
        y_vals[k+2] = sum(derivative(i) * ((x_k - x0)**(i)) /
                          math.factorial(i) for i in range(6))

    return x_vals, y_vals


# Метод Рунге-Кутты 4-го порядка
def runge_kutta_4(x0, y0, h, N):
    x_vals = [x0]
    y_vals = [y0]
    for k in range(N):
        x_k = x_vals[-1]
        y_k = y_vals[-1]
        k1 = h * f(x_k, y_k)
        k2 = h * f(x_k + h / 2, y_k + k1 / 2)
        k3 = h * f(x_k + h / 2, y_k + k2 / 2)
        k4 = h * f(x_k + h, y_k + k3)
        y_next = y_k + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x_vals.append(x_k + h)
        y_vals.append(y_next)
    return x_vals, y_vals


# Метод Эйлера
def euler_method(x0, y0, h, N):
    x_vals = [x0]
    y_vals = [y0]
    for k in range(N):
        x_k = x_vals[-1]
        y_k = y_vals[-1]
        y_next = y_k + h * f(x_k, y_k)
        x_vals.append(x_k + h)
        y_vals.append(y_next)
    return x_vals, y_vals


# Метод Эйлера I (средний прямоугольник)
def euler_method_1(x0, y0, h, N):
    x_vals = [x0]
    y_vals = [y0]
    for k in range(N):
        x_k = x_vals[-1]
        y_k = y_vals[-1]
        y_half = y_k + (h / 2) * f(x_k, y_k)
        y_next = y_k + h * f(x_k + h / 2, y_half)
        x_vals.append(x_k + h)
        y_vals.append(y_next)
    return x_vals, y_vals


# Метод Эйлера II (трапеция)
def euler_method_2(x0, y0, h, N):
    x_vals = [x0]
    y_vals = [y0]
    for k in range(N):
        x_k = x_vals[-1]
        y_k = y_vals[-1]
        y_predict = y_k + h * f(x_k, y_k)
        y_next = y_k + (h / 2) * (f(x_k, y_k) + f(x_k + h, y_predict))
        x_vals.append(x_k + h)
        y_vals.append(y_next)
    return x_vals, y_vals


def adams_4th_order(x0, y0, h, N):
    # Получаем начальные значения из метода Тейлора
    x, y = taylor_method(x0, y0, h, N)

    for i in range(4, len(y)-1):
        y[i+1] = y[i] + h/24 * (55 * f(x[i], y[i]) - 59 * f(x[i-1],
                                y[i-1]) + 37 * f(x[i-2], y[i-2]) - 9 * f(x[i-3], y[i-3]))
        # print(y[i+1])
    return x, y


def adjust_table(arr, max_len):
    while (len(arr) < max_len):
        arr.insert(0, None)


def main():
    while True:
        try:
            print("\nВведите параметры задачи:")
            h = float(input("Введите шаг h: "))
            N = int(input("Введите количество шагов N: "))
            x0, y0 = 0, 0.25

            x_t, y_taylor = taylor_method(x0, y0, h, N)
            x_rk, y_rk4 = runge_kutta_4(x0, y0, h, N)
            x_e, y_euler = euler_method(x0, y0, h, N)
            x_e1, y_euler1 = euler_method_1(x0, y0, h, N)
            x_e2, y_euler2 = euler_method_2(x0, y0, h, N)
            x_adams, y_adams = adams_4th_order(x0, y0, h, N)

            max_len = max(len(y_taylor), len(y_rk4), len(y_euler),
                          len(y_euler1), len(y_euler2), len(y_adams))
            adjust_table(y_taylor, max_len)
            adjust_table(y_rk4, max_len)
            adjust_table(y_euler, max_len)
            adjust_table(y_euler1, max_len)
            adjust_table(y_euler2, max_len)

            table_data = []
            for i in range(max_len):
                exact = exact_solution(x_t[i])
                table_data.append([
                    i-2, x_t[i], exact,
                    y_taylor[i], str(
                        f'{abs(exact-y_taylor[i]):.13f}'), y_rk4[i], y_euler[i],
                    y_euler1[i], y_euler2[i], y_adams[i]
                ])

            last_idx = max_len - 1
            exact = exact_solution(x_t[last_idx])
            table_data.append([
                "Погрешности",
                None, None,
                str(f'{abs(y_taylor[last_idx] - exact):.13f}'), None,
                str(f'{abs(y_rk4[last_idx] - exact):.13f}'),
                str(f'{abs(y_euler[last_idx] - exact):.13f}'),
                str(f'{abs(y_euler1[last_idx] - exact):.13f}'),
                str(f'{abs(y_euler2[last_idx] - exact):.13f}'),
                str(f'{abs(y_adams[last_idx] - exact):.13f}')
            ])

            headers = [
                "Шаг", "x", "Точное",
                "Тейлор", "Погрешность Тэйлор", "Рунге-Кутта", "Эйлер", "Эйлер I", "Эйлер II", "Адамс 4-го порядка"
            ]

            print("\nТаблица решений:")
            print(tabulate(table_data, headers=headers,
                  floatfmt=".11f", tablefmt="grid"))

            choice = input(
                "\nВыберите действие:\n1. Ввести новые параметры N и h.\n2. Завершить программу.\nВаш выбор: ")
            if choice == "1":
                continue
            elif choice == "2":
                break
            else:
                print("Неверный выбор, попробуйте снова.")
        except ValueError as e:
            print(f"Ошибка ввода: {e}. Попробуйте снова.")


if __name__ == "__main__":
    main()
