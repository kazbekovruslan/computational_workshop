import math


def f(x): return x * math.sin(x) - 1


def df(x): return math.sin(x) + x * math.cos(x)


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
    while b - a > 2 * epsilon:
        steps += 1
        c = (a + b) / 2
        if (f(a) * f(c) <= 0):
            b = c
        else:
            a = c
    return (a + b) / 2, steps


def newton_method(f, df, x0, eps):
    steps = 0
    while True:
        x1 = x0 - f(x0) / df(x0)
        if abs(x1 - x0) < eps:
            return x1, steps
        x0 = x1
        steps += 1


def modified_newton_method(f, df, x0, eps):
    steps = 0
    df_x0 = df(x0)
    while True:
        x1 = x0 - f(x0) / df_x0
        if abs(x1 - x0) < eps:
            return x1, steps
        x0 = x1
        steps += 1


def secant_method(f, a, b, eps):
    steps = 0
    while abs(b - a) >= eps:
        x2 = b - f(b) * (b - a) / (f(b) - f(a))
        a, b = b, x2
        steps += 1
    return b, steps


def main():
    print("ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ")
    print(
        "Данные в условии: f(x)= x * sin(x) - 1 [A, B] = [-10; 2] ε = 10^-5")

    A, B = int(input("Введите параметр A: ")), int(
        input("Введите параметр B: "))
    epsilon = 10 ** ((-1) *
                     int(input("Введите x -- степень точности 10 ** (-x): ")))
    N = 100

    intervals = find_intervals(A, B, f, N)

    results = []

    for a, b in intervals:
        bisection_result, b_steps = bisection_method(f, a, b, epsilon)
        newton_result, n_steps = newton_method(f, df, a, epsilon)
        modified_newton_result, n_steps = modified_newton_method(
            f, df, a, epsilon)
        secant_result, s_steps = secant_method(f, a, b, epsilon)
        results.append({
            "interval": (a, b),
            "bisection": (bisection_result, b_steps),
            "newton": (newton_result, n_steps),
            "modified_newton": (modified_newton_result, n_steps),
            "secant": (secant_result, s_steps)
        })

    for result in results:
        print("\nИнтервал: ", result['interval'])
        print("Метод бисекции: ", result['bisection'])
        print("Метод Ньютона: ", result['newton'])
        print("Модифицированный метод Ньютона: ",
              result['modified_newton'])
        print("Метод секущих: ", result['secant'])


if __name__ == "__main__":
    main()
