import numpy as np
import math
from tabulate import tabulate
from scipy.integrate import quad


def f(x):
    return math.sin(x)


def exact_integral(a, b):
    def F(x):
        return -0.5 * math.exp(-x) * (math.sin(x) + math.cos(x))

    return F(b) - F(a)


"""[a,b] ∫ x^k *\rho(x) dx"""


def weight_moments(a, b, degree):
    moments = []
    for k in range(degree + 1):
        moment, _ = quad(lambda x: x**k * math.exp(-x), a, b)
        moments.append(moment)
    return moments


"""чтобы построить ИКФ, она должна быть точной для всех многочленов степени N-1. Это приводит к системе линейных уравнений:
∑_{k=1}^N Aₖ xₖʲ = μⱼ,  j = 0, 1, …, N-1.

где:
    μⱼ = ∫ₐᵇ xʲ * rho(x) dx, где j = 0, 1, 2, ..., N-1.
 — j-й момент весовой функции."""


def solve_coefficients(nodes, moments):

    n = len(nodes)
    A = np.vander(nodes, increasing=True).T[:n]
    b = moments[:n]
    coefficients = np.linalg.solve(A, b)
    return coefficients


def calculate_ikf_integral(f_values, coefficients):
    return np.dot(coefficients, f_values)  # скалярное произведение


def check_accuracy(coefficients, nodes, N, a, b):
    for degree in range(1, N):  # Проверяем для многочленов степени 1, 2, ..., N-1

        def poly(x, degree=degree):
            return x**degree

        # Вычисляем точное значение интеграла для многочлена
        exact_value, _ = quad(lambda x: poly(x, degree) * math.exp(-x), a, b)

        # Вычисляем приближенное значение через ИКФ
        f_values = [poly(x) for x in nodes]
        ikf_value = calculate_ikf_integral(f_values, coefficients)

        # Сравниваем точное и приближенное значение
        print(f"Для многочлена степени {degree}:")
        print(f"Точное значение интеграла: {exact_value:.12f}")
        print(f"Приближенное значение ИКФ: {ikf_value:.12f}")
        print(f"Ошибка: {abs(exact_value - ikf_value):.12e}\n")


def main():
    print(
        "Задача 4.1: Приближённое вычисление интегралов при помощи интерполяционных квадратурных формул (ИКФ)"
    )
    print("Вариант 13: f(x) = sin(x), rho(x)=e^(-x)")
    print("Введите границы интегрирования")
    a, b = map(float, input().split())

    print(f"Отрезок интегрирования: [{a}, {b}]")

    exact_value = exact_integral(a, b)
    print(f"Точное значение интеграла: {exact_value:.12f}")

    N = int(input("Введите количество узлов (N): "))
    print("Введите узлы (через пробел):")

    nodes = list(map(float, input().split()))

    def check_in(nodes):
        check = False
        for i in nodes:
            check = check or (i < a or i > b)
        return check

    check = check_in(nodes)

    while len(set(nodes)) != N or check:
        print(
            f"Узлы должны быть попарно различны, лежать в [{a},{b}] и их должно быть {N}"
        )
        print(f"Введите {N} узлов (через пробел):")
        nodes = list(map(float, input().split()))
        check = check_in(nodes)

    # Шаг 3: Вычисление моментов весовой функции
    moments = weight_moments(a, b, N - 1)
    print("\nМоменты весовой функции:")
    for i, moment in enumerate(moments):
        print(f"μ_{i} =  ∫ₐᵇ x^{i} * e^(-x) dx = {moment:.12f}")

    # Шаг 4: Решение СЛАУ для коэффициентов
    coefficients = solve_coefficients(nodes, moments)
    print("\nТаблица узлов и коэффициентов:")
    table = [[f"x_{i+1}", nodes[i], coefficients[i]] for i in range(N)]
    print(
        tabulate(
            table, headers=["Узел", "Значение узла", "Коэффициент"], floatfmt=".12f"
        )
    )
    print("Сумма коэффициентов:", sum(coefficients))
    print()
    print("Проверка точности:")
    check_accuracy(coefficients, nodes, N, a, b)
    print("Результат решения:")
    f_values = [f(node) for node in nodes]
    ikf_value = calculate_ikf_integral(f_values, coefficients)
    print(f"Значение интеграла с помощью ИКФ: {ikf_value:.12f}")
    print(f"Точное значение: {exact_value:.12f}")
    print(f"Разница между точным и ИКФ: {abs(exact_value - ikf_value):.12e}")


if __name__ == "__main__":
    main()
