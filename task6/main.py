from math import sin, sqrt, exp
from scipy import integrate
from scipy.linalg import solve
import numpy as np


def ro(x):
    return pow(x, -0.25)


def gauss_ro(x):
    return 1


def meler_ro(x):
    return 1 / sqrt(1 - x * x)


def f(x):
    return sin(x)


def get_input(prompt, type_fn=float):
    return type_fn(input(prompt))


def calculate_moments(weight_function, a, b, N):
    moments = []
    for deg in range(N):
        moment = integrate.quad(lambda x: x ** deg *
                                weight_function(x), a, b)[0]
        print(f' m_{deg} = {moment:.12f}')
        moments.append(moment)
    return moments


def calculate_coefs(points, moments):
    N = len(points)
    LHS = [[x ** deg for x in points] for deg in range(N)]
    return solve(np.array(LHS), np.array(moments).reshape(N, 1))


def calculate_integral(coefs, points, func):
    return sum(coef * func(point) for coef, point in zip(coefs.flatten(), points))


def main():
    print('Приближённое вычисление интегралов при помощи квадратурных формул Наивысшей Алгебраической Степени Точности (КФ НАСТ)')
    print('∫p(x)f(x). p=x^-0.25, f=sin(x)')
    a_const = 0
    b_const = 1

    a, b = get_input("Введите A: "), get_input("Введите B: ")
    N = get_input("Введите количество узлов: ", int)

    precise_integral, _ = integrate.quad(lambda x: ro(x) * f(x), a, b)

    points = np.linspace(a, b, N)
    # Форматирование узлов для вывода
    moments = calculate_moments(ro, a, b, N)

    print('\nУзлы и их коэффициенты:')
    coefs = calculate_coefs(points, moments)
    for i, point in enumerate(points):
        print(f'Узел {point:.12f}: коэффициент = {coefs[i, 0]:.12f}')

    def test_polynom(x): return x ** (N - 1) + 1
    lhs = integrate.quad(lambda x: ro(x) * test_polynom(x), a, b)[0]
    rhs = calculate_integral(coefs, points, test_polynom)
    print('\nПроверка точности:')
    print(f'Многочлен x^{N - 1}: {abs(lhs - rhs):.12e}')

    iqf_result = calculate_integral(coefs, points, f)
    print(f'\nТочное значение: {precise_integral:.12f}')
    print(f'Значение интеграла с помощью ИКФ: {iqf_result:.12f}')
    print(
        f'Разница между точным и ИКФ: {abs(precise_integral - iqf_result):.12e}')

    print('=== КФНАСТ ===')

    weight_functions = [ro, gauss_ro, meler_ro]

    sections = [(a_const, b_const), (-1, 1), (-1, 1)]

    ro_strings = ['x^-0.25', '1', '1/sqrt(1 - x^2)']
    w_polys = ['x^2-1.06682x+0.20284', '1.5 * x^2 - 0.5', '2x^2 - 1']
    roots_of_ws = [[0.247601973345, 0.819218026654],
                   [-1 / sqrt(3), 1 / sqrt(3)],
                   [-1 / sqrt(2), 1 / sqrt(2)]]
    degree_of_w = 2

    for i, (weight_function, (left, right), ro_string, w_poly, roots_of_w) in enumerate(
            zip(weight_functions, sections, ro_strings, w_polys, roots_of_ws)):
        precise_integral, _ = integrate.quad(
            lambda x: weight_function(x) * f(x), left, right)
        print(f'\np = {ro_string}:\n')

        print('Moments of weight function:')
        moments = calculate_moments(weight_function, left, right, degree_of_w)

        print(f'Orthogonal poly: {w_poly}')
        coefs = calculate_coefs(roots_of_w, moments)
        print('\nNodes and corresponding coefficients:')
        for i, root in enumerate(roots_of_w):
            print(f'  Node {root:.12f}: coeff = {coefs[i, 0]:.12f}')

        def test_polynom(x): return x ** 3 + x ** 2 - x - 1
        lhs = integrate.quad(lambda x: weight_function(
            x) * test_polynom(x), left, right)[0]
        rhs = calculate_integral(coefs, roots_of_w, test_polynom)
        print('\n Accuracy check on polynomials:')
        print(f'Polynom x^3 + x^2 - x - 1: {abs(lhs - rhs):.12e}')

        iqf_result = calculate_integral(coefs, roots_of_w, f)
        print(f'\nValue КФ НАСТ: {iqf_result:.12f}')
        print(f'Exact: {precise_integral:.12f}')
        print(f'Difference: {abs(precise_integral - iqf_result):.12e}')


if __name__ == "__main__":
    main()
