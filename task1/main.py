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
    return (a + b) / 2, steps, first_approximation, b - a


def find_first_approximation(a, b, f, ddf):
    x0 = random.uniform(a, b)
    while f(x0) * ddf(x0) <= 0:
        x0 = random.uniform(a, b)
    return x0


# change
def newton_method(f, df, ddf, a, b, epsilon):
    steps = 0
    first_approximation = find_first_approximation(a, b, f, ddf)
    x0 = first_approximation
    x1 = x0 - f(x0) / df(x0)
    while abs(x1 - x0) >= epsilon:
        x1 = x0 - f(x0) / df(x0)
        x0 = x1
        steps += 1
    return x1, steps, first_approximation, abs(x1 - x0)


def modified_newton_method(f, df, ddf, a, b, epsilon):
    steps = 0
    first_approximation = find_first_approximation(a, b, f, ddf)
    x0 = first_approximation
    df_x0 = df(x0)
    x1 = x0 - f(x0) / df_x0
    while abs(x1 - x0) >= epsilon:
        x1 = x0 - f(x0) / df_x0
        x0 = x1
        steps += 1
    return x1, steps, first_approximation, abs(x1 - x0)


def secant_method(f, a, b, eps):
    steps = 0
    while abs(b - a) >= eps:
        x2 = b - f(b) * (b - a) / (f(b) - f(a))
        a, b = b, x2
        steps += 1
    return b, steps, a, abs(b - a)


def main():
    print("=== ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ ===")
    print("Уравнение: f(x) = x * sin(x) - 1")
    print("Интервал: [A, B] = [-10, 2]")
    print("Точность: ε = 10^-x")

    A = int(input("Введите параметр A: "))
    B = int(input("Введите параметр B: "))
    epsilon = 10 ** (-int(input("Введите x (степень точности): ")))
    N = 100

    print("\n=== Поиск интервалов изменения знака ===")
    intervals = find_intervals(A, B, f, N)
    print(f"Найдено {len(intervals)} интервалов: {intervals}")

    print("\n=== Уточнение корней ===")
    for idx, (a, b) in enumerate(intervals, start=1):
        print(f"\nИнтервал {idx}: [{a:.6f}, {b:.6f}]")
        bisection_result, b_steps, b_first_approximation, b_diff = bisection_method(
            f, a, b, epsilon)
        newton_result, n_steps, n_first_approximation, n_diff = newton_method(
            f, df, ddf, a, b, epsilon)
        modified_newton_result, mn_steps, mn_first_approximation, mn_diff = modified_newton_method(
            f, df, ddf, a, b, epsilon)
        secant_result, s_steps, s_first_approximation, s_diff = secant_method(
            f, a, b, epsilon)

        print(f"Метод бисекции:")
        print(
            f"  Корень: {bisection_result:.6f}, Шагов: {b_steps}, Начальное приближение: {b_first_approximation:.6f}, Расстояние до предыдущего приближения: {abs(b_diff)}")
        print(f"Метод Ньютона:")
        print(
            f"  Корень: {newton_result:.6f}, Шагов: {n_steps}, Начальное приближение: {n_first_approximation:.6f}, Расстояние до предыдущего приближения: {abs(n_diff)}")
        print(f"Модифицированный метод Ньютона:")
        print(
            f"  Корень: {modified_newton_result:.6f}, Шагов: {mn_steps}, Начальное приближение: {mn_first_approximation:.6f}, Расстояние до предыдущего приближения: {abs(mn_diff)}")
        print(f"Метод секущих:")
        print(
            f"  Корень: {secant_result:.6f}, Шагов: {s_steps}, Начальное приближение: {s_first_approximation:.6f}, Расстояние до предыдущего приближения: {abs(s_diff)}")

    print("\n=== Завершено! ===")


if __name__ == "__main__":
    main()
