import numpy as np
import math
from tabulate import tabulate
from scipy.integrate import quad
from scipy.linalg import solve


# сумма A_k -- нулевой момент функции


def f(x):
    return math.sin(x)


def print_table(table):
    if not table or not table[0]:
        print("Таблица пуста или имеет неверный формат.")
        return

    # Рассчитываем ширину каждой колонки
    col_widths = [max(len(str(row[i])) for row in table)
                  for i in range(len(table[0]))]

    # Печать строки разделителя
    def print_separator():
        print("+" + "+".join("-" * (width + 2) for width in col_widths) + "+")

    # Печать строки с данными (альтернативный способ)
    def print_row(row):
        row_str = "| "
        for i, cell in enumerate(row):
            row_str += f"{str(cell):<{col_widths[i]}} | "
        print(row_str)

    # Печатаем таблицу
    print_separator()  # Верхний разделитель
    print_row(table[0])  # Заголовок таблицы
    print_separator()  # Разделитель под заголовком
    for row in table[1:]:
        print_row(row)  # Строки с данными
    print_separator()  # Нижний разделитель


def exact_integral(a, b):
    def F(x):
        return 0.5 * math.exp(x) * (math.sin(x) - math.cos(x))

    return F(b) - F(a)


def weight_moments(a, b, degree):
    moments = []
    for k in range(degree + 1):
        moment, _ = quad(lambda x: x**k * math.exp(x), a, b)
        moments.append(moment)
    return moments


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
        exact_value, _ = quad(lambda x: poly(x, degree) * math.exp(x), a, b)

        # Вычисляем приближенное значение через ИКФ
        f_values = [poly(x) for x in nodes]
        ikf_value = calculate_ikf_integral(f_values, coefficients)

        # Сравниваем точное и приближенное значение
        print(f"Многочлен степени {degree}:")
        print(f"Точное значение интеграла: {exact_value:.12f}")
        print(f"Приближенное значение ИКФ: {ikf_value:.12f}")
        print(f"Погрешность: {abs(exact_value - ikf_value):.12e}\n")


def get_input(prompt, type_fn=float):
    return type_fn(input(prompt))


def calculate_moments(weight_function, a, b, N):
    moments = []
    for deg in range(N):
        moment = quad(lambda x: x ** deg *
                      weight_function(x), a, b)[0]
        # print(f' m_{deg} = {moment:.12f}')
        moments.append(moment)
    return moments


def calculate_coefs(points, moments):
    N = len(points)
    LHS = [[x ** deg for x in points] for deg in range(N)]
    return solve(np.array(LHS), np.array(moments).reshape(N, 1))


def calculate_integral(coefs, points, func):
    return sum(coef * func(point) for coef, point in zip(coefs.flatten(), points))


def main():
    print(
        "Задача 6: Приближённое вычисление интегралов при помощи квадратурных формул Наивысшей Алгебраической Степени Точности (КФ НАСТ)"
    )
    print("f(x) = sin(x), ρ(x)=e^x")
    a = float(input("Введите левую границу интегрирования: "))
    b = float(input("Введите правую границу интегрирования: "))

    print(f"Отрезок интегрирования: [{a}, {b}]")

    exact_value = exact_integral(a, b)
    print(f"Точное значение интеграла: {exact_value:.12f}")

    N = int(input("\nВведите количество узлов ИКФ (N): "))
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
        print(f"μ_{i} =  [a,b] ∫ x^{i} * e^x dx = {moment:.12f}")

    # Шаг 4: Решение СЛАУ для коэффициентов
    coefficients = solve_coefficients(nodes, moments)

    print("\nТаблица узлов и коэффициентов:")
    table = [["№", "Узел", "Коэффициент"]]
    table_1 = [[f"x_{i+1}", nodes[i], coefficients[i]] for i in range(N)]
    table += table_1

    print_table(table)

    print("Сумма коэффициентов:", sum(coefficients))
    print()
    print("Проверка точности:")
    check_accuracy(coefficients, nodes, N, a, b)

    f_values = [f(node) for node in nodes]
    ikf_value = calculate_ikf_integral(f_values, coefficients)
    print(f"\nЗначение интеграла с помощью ИКФ: {ikf_value:.12f}")
    print(f"Точное значение: {exact_value:.12f}")
    print(f"Разница между точным и ИКФ: {abs(exact_value - ikf_value):.12e}")

    print('=== КФНАСТ ===')

    def ro(x):
        return math.exp(x)

    def gauss_ro(x):
        return 1

    def meler_ro(x):
        return 1 / math.sqrt(1 - x * x)

    a_const = 0
    b_const = 1

    weight_functions = [ro, gauss_ro, meler_ro]

    sections = [(a_const, b_const), (-1, 1), (-1, 1)]

    ro_strings = ['e^x', '1', '1/sqrt(1 - x^2)']
    w_polys = ['x^2-1.06682x+0.20284', '1.5 * x^2 - 0.5', '2x^2 - 1']
    roots_of_ws = [[0.247601973345, 0.819218026654],
                   [-1 / math.sqrt(3), 1 / math.sqrt(3)],
                   [-1 / math.sqrt(2), 1 / math.sqrt(2)]]
    degree_of_w = 2

    for i, (weight_function, (left, right), ro_string, w_poly, roots_of_w) in enumerate(
            zip(weight_functions, sections, ro_strings, w_polys, roots_of_ws)):
        precise_integral, _ = quad(
            lambda x: weight_function(x) * f(x), left, right)
        print(f'\n=== p = {ro_string} ===\n')

        # print('Моменты весовой функции:')
        moments = calculate_moments(weight_function, left, right, degree_of_w)

        print(f'Ортогональный многочлен: {w_poly}')
        coefs = calculate_coefs(roots_of_w, moments)
        # print('\nУзлы и коэффициенты:')
        # for i, root in enumerate(roots_of_w):
        #     print(f'  Узел {root:.12f}: коэффициент = {coefs[i, 0]:.12f}')

        print("\nМоменты весовой функции:")
        for i, moment in enumerate(moments):
            # print(f"μ_{i} =  [a,b] ∫ x^{i} * e^x dx = {moment:.12f}")
            print(f"μ_{i} = {moment:.12f}")

        print("\nТаблица узлов и коэффициентов:")
        table = [["№", "Узел", "Коэффициент"]]
        table_1 = [[f"x_{i+1}", roots_of_w[i], coefs[i, 0]]
                   for i in range(degree_of_w)]
        table += table_1

        print_table(table)

        def test_polynom(x): return x ** 3 + x ** 2 - x - 1
        lhs = quad(lambda x: weight_function(
            x) * test_polynom(x), left, right)[0]
        rhs = calculate_integral(coefs, roots_of_w, test_polynom)
        print('\n Проверка точности:')
        print(
            f'Многочлен x^3 + x^2 - x - 1\nПогрешность: {abs(lhs - rhs):.12e}')

        iqf_result = calculate_integral(coefs, roots_of_w, f)
        print(f'\nЗначение КФ НАСТ: {iqf_result:.12f}')
        print(f'Точное значение: {precise_integral:.12f}')
        print(f'Разница: {abs(precise_integral - iqf_result):.12e}')


if __name__ == "__main__":
    main()
