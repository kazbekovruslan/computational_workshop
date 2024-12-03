import math
import random


def f(x): return x * math.sin(x) - 1


def df(x): return math.sin(x) + x * math.cos(x)


def ddf(x): return 2 * math.cos(x) - x * math.sin(x)


def find_intervals(a, b, f, n):
    intervals = []
    h = (b - a) / n

    begin, end = a, a + h
    while end <= b:
        if f(begin) * f(end) <= 0:
            intervals.append((begin, end))
        begin, end = end, end + h
    return intervals


def bisection_method(f, a, b, epsilon):
    steps = 0
    first_approximation = (a + b) / 2
    while b - a > 2 * epsilon:
        steps += 1
        c = (a + b) / 2
        if (f(a) * f(c) <= 0):
            b = c
        else:
            a = c
    xm = (a + b) / 2
    return xm, steps, first_approximation, b - a, abs(f(xm))


def find_first_approximation(a, b, f, ddf):
    x0 = random.uniform(a, b)
    while f(x0) * ddf(x0) <= 0:
        x0 = random.uniform(a, b)
    return x0


def newton_method(f, df, ddf, a, b, epsilon):
    steps = 0
    first_approximation = find_first_approximation(a, b, f, ddf)
    x0 = first_approximation
    x1 = x0 - f(x0) / df(x0)
    while abs(x1 - x0) >= epsilon:
        x1, x0 = x0, x0 - f(x0) / df(x0)
        steps += 1
    return x1, steps, first_approximation, abs(x1 - x0), abs(f(x1))


def modified_newton_method(f, df, ddf, a, b, epsilon):
    steps = 0
    first_approximation = find_first_approximation(a, b, f, ddf)
    x0 = first_approximation
    df_x0 = df(x0)
    x1 = x0 - f(x0) / df_x0
    while abs(x1 - x0) >= epsilon:
        x1, x0 = x0, x0 - f(x0) / df_x0
        steps += 1
    return x1, steps, first_approximation, abs(x1 - x0), abs(f(x1))


def secant_method(f, a, b, eps):
    steps = 0
    first_approximation = a
    while abs(b - a) >= eps:
        x2 = b - f(b) * (b - a) / (f(b) - f(a))
        a, b = b, x2
        steps += 1
    return b, steps, first_approximation, abs(b - a), abs(f(b))


def main():
    print("=== ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ ===")
    print("Вариант №5")
    print("Функция: f(x) = x * sin(x) - 1")
    print("Интервал: [A, B] = [-10, 2]")
    print("Точность: ε = 10^{-5}")

    A = int(input("Введите параметр A: "))
    B = int(input("Введите параметр B: "))
    epsilon = 10 ** (-int(input("Введите степень точности ε: ")))
    N = 100

    otd = True
    while otd:
        N = int(input("Введите на сколько отрезков поделится изначальный отрезок: "))
        print(f"Длина интервала = {(B - A) / N}")

        print("\n=== Поиск интервалов изменения знака ===")
        intervals = find_intervals(A, B, f, N)
        print(f"Найдено {len(intervals)} интервалов:")

        for i, (start, end) in enumerate(intervals, start=1):
            print(f"  Интервал {i}: [{start:.3f}, {end:.3f}]")

        answer = input(
            "Хотите изменить количество отрезков разбиения? (yes or no): ")
        if answer == "no":
            otd = False

    print("\n=== Уточнение корней ===")
    results = []

    for idx, (a, b) in enumerate(intervals, start=1):
        print(f"\n=== Интервал {idx}: [{a:.3f}, {b:.3f}] ===")

        b_result, b_steps, b_first_approximation, b_diff, b_residual = bisection_method(
            f, a, b, epsilon)

        n_result, n_steps, n_first_approximation, n_diff, n_residual = newton_method(
            f, df, ddf, a, b, epsilon)

        mn_result, mn_steps, mn_first_approximation, mn_diff, mn_residual = modified_newton_method(
            f, df, ddf, a, b, epsilon)

        s_result, s_steps, s_first_approximation, s_diff, s_residual = secant_method(
            f, a, b, epsilon)

        methods = [
            ("Метод бисекции", b_first_approximation,
                b_steps, b_result, b_diff, b_residual),
            ("Метод Ньютона", n_first_approximation,
                n_steps, n_result, n_diff, n_residual),
            ("Модифицированный метод Ньютона", mn_first_approximation,
                mn_steps, mn_result, mn_diff, mn_residual),
            ("Метод секущих", s_first_approximation,
                s_steps, s_result, s_diff, s_residual)
        ]

        for method_name, first_approximation, steps, result, diff, residual in methods:
            print(f"{method_name}:")
            print(f"  Начальное приближение: {first_approximation:}")
            print(f"  Количество шагов: {steps}")
            print(f"  Найденное решение: {result:}")
            print(f"  |x_m - x_m-1|: {diff}")
            print(f"  Невязка: {residual}")

        results.append((idx, methods))
    print("=== Расчет завершен ===")


if __name__ == "__main__":
    main()
