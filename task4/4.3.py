import numpy as np
import math


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


def f(x):
    return math.sin(x)


def f_0(x):
    return 1


def f_1(x):
    return 2 * x


def f_2(x):
    return 3 * (x**2)


def f_3(x):
    return 4 * (x**3)


func_names = {
    f: "sin(x)",
    f_0: "1",
    f_1: "2x",
    f_2: "3x^2",
    f_3: "4x^3",
}


def func_to_str(func):
    return func_names.get(func, "Unknown function")


def exact_integral(f, a, b):
    from scipy.integrate import quad

    result, _ = quad(f, a, b)
    return result


def composite_quadratures(f, a, b, m):
    def complex_left_rectangle_method(f, a, b, m):
        h = (b - a) / m
        return h * sum(f(a + i * h) for i in range(0, m))

    def complex_right_rectangle_method(f, a, b, m):
        h = (b - a) / m
        return h * sum(f(a + i * h) for i in range(1, m + 1))

    def complex_middle_rectangle_method(f, a, b, m):
        h = (b - a) / m
        return h * sum(f(a + i * h + h / 2) for i in range(0, m))

    def complex_trapeze_method(f, a, b, m):
        return (
            complex_left_rectangle_method(f, a, b, m)
            + complex_right_rectangle_method(f, a, b, m)
        ) / 2

    def complex_simpson_method(f, a, b, m):
        return (
            complex_left_rectangle_method(f, a, b, m)
            + complex_right_rectangle_method(f, a, b, m)
            + 4 * complex_middle_rectangle_method(f, a, b, m)
        ) / 6

    left_rect = complex_left_rectangle_method(f, a, b, m)
    right_rect = complex_right_rectangle_method(f, a, b, m)
    mid_rect = complex_middle_rectangle_method(f, a, b, m)
    trapezoid = complex_trapeze_method(f, a, b, m)
    simpson = complex_simpson_method(f, a, b, m)

    if m == 1:
        x = np.linspace(a, b, m + 1)
        three_eighths = 0
        for i in range(m):
            a_k = x[i]
            b_k = x[i + 1]
            h_38 = (b_k - a_k) / 3
            three_eighths += (b_k - a_k) * (
                (1 / 8) * f(a_k) +
                (3 / 8) * f(a_k + h_38) +
                (3 / 8) * f(a_k + 2 * h_38) +
                (1 / 8) * f(b_k)
            )
        return left_rect, right_rect, mid_rect, trapezoid, simpson, three_eighths

    return left_rect, right_rect, mid_rect, trapezoid, simpson


def runge_correction_pair(j_h, j_h2, p):
    return j_h2 + (j_h2 - j_h) / (2**p - 1)


def main():
    print("Приближённое вычисление интеграла по составным квадратурным формулам.")

    running = True
    while running:
        function = int(
            input(
                "Выберите функцию:\n1. f = sin(x)\n2. f = 1\n3. f = 2x\n4. f = 3x^2\n5. f = 4x^3\n"
            )
        )
        while function not in [1, 2, 3, 4, 5]:
            print("Некорректный ввод! Выберите функцию (введите цифру от 1 до 4)")
            function = int(
                input(
                    "Выберите функцию:\n1. f = sin(x)\n2. f = 1\n3. f = 2x\n4. f = 3x^2\n5. f = 4x^3\n"
                )
            )

        match function:
            case 1:
                func = f
            case 2:
                func = f_0
            case 3:
                func = f_1
            case 4:
                func = f_2
            case 5:
                func = f_3

        a = float(input("Введите нижний предел интегрирования A: "))
        b = float(input("Введите верхний предел интегрирования B: "))
        m = int(input("Введите число промежутков деления m: "))

        j_exact = exact_integral(func, a, b)
        h = (b - a) / m
        print(f"\nТочное значение интеграла: J = {j_exact:.15f}")
        print(
            f"Параметры задачи: A={a}, B={b}, m={m}, h={h:.15f}, f = {func_to_str(func)}\n"
        )

        j_approx = composite_quadratures(func, a, b, m)
        j_approx_h2 = composite_quadratures(func, a, b, m * 2)

        methods = [
            "Левые прямоугольники",
            "Правые прямоугольники",
            "Средние прямоугольники",
            "Трапеции",
            "Симпсона",
        ]
        orders = [1, 1, 2, 2, 4]

        if m == 1:
            methods = [
                "Левые прямоугольники",
                "Правые прямоугольники",
                "Средние прямоугольники",
                "Трапеции",
                "Метод 3/8",
                "Симпсона",
            ]

        table = [["Метод", "J(h)", "Абс. погр.", "Отн. погр."]]
        for method, j in zip(methods, j_approx):
            abs_error = abs(j_exact - j)
            rel_error = abs_error / abs(j_exact)
            table.append(
                [method, f"{j:.15f}", f"{abs_error:.15e}", f"{rel_error:.15e}"]
            )

        print(f"Результаты для m = {m}:")
        print_table(table)

        if m == 1:
            print("=== Завершение работы ===")
            return
        # Таблица результатов с уточнением по Рунге для m
        table_runge = [["Метод", "J(точн.)", "Абс. погр.", "Отн. погр."]]
        for method, j, j_h2, p in zip(methods, j_approx, j_approx_h2, orders):
            j_refined = runge_correction_pair(j, j_h2, p)
            abs_error_ref = abs(j_exact - j_refined)
            rel_error_ref = abs_error_ref / abs(j_exact)
            table_runge.append(
                [
                    method,
                    f"{j_refined:.15f}",
                    f"{abs_error_ref:.15e}",
                    f"{rel_error_ref:.15e}",
                ]
            )

        print(f"\nУточнённые значения методом Рунге для m = {m}:")
        print_table(table_runge)

        # Уточнение для m*l
        l = int(input("\nВведите параметр l для увеличения числа разбиений: "))
        j_approx_l = composite_quadratures(func, a, b, m * l)
        j_approx_h2l = composite_quadratures(func, a, b, m * 2 * l)

        # Таблица результатов без уточнения для m*l
        table_l = [["Метод", "J(h/l)", "Абс. погр.", "Отн. погр."]]
        for method, j_l in zip(methods, j_approx_l):
            abs_error_l = abs(j_exact - j_l)
            rel_error_l = abs_error_l / abs(j_exact)
            table_l.append(
                [method, f"{j_l:.15f}", f"{abs_error_l:.15e}",
                    f"{rel_error_l:.15e}"]
            )

        print(f"\nРезультаты для m*l= {m*l}:")
        print_table(table_l)

        # Таблица результатов с уточнением по Рунге для m*l
        table_runge_l = [
            ["Метод", "J(точн. для m*l)", "Абс. погр.", "Отн. погр."]]
        for method, j_l, j_h2l, p in zip(methods, j_approx_l, j_approx_h2l, orders):
            j_specified_l = runge_correction_pair(j_l, j_h2l, p)
            abs_error_ref_l = abs(j_exact - j_specified_l)
            rel_error_ref_l = abs_error_ref_l / abs(j_exact)
            table_runge_l.append(
                [
                    method,
                    f"{j_specified_l:.15f}",
                    f"{abs_error_ref_l:.15e}",
                    f"{rel_error_ref_l:.15e}",
                ]
            )

        print(f"\nУточнённые значения методом Рунге для m*l = {m*l}:")
        print_table(table_runge_l)

        answer = int(
            input(
                "Введите 1, чтобы повторить с другими параметрами; введите 2, чтобы выйти: "
            )
        )
        while answer not in [1, 2]:
            print("Некорректный ввод! Введите 1 или 2")
            answer = int(
                input(
                    "Введите 1, чтобы повторить с другими параметрами; введите 2, чтобы выйти: "
                )
            )
        if answer == 2:
            running = False
            print("=== Работа программы завершена ===")


if __name__ == "__main__":
    main()
