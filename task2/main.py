import math


def f(x):
    return 1 - math.exp(-2 * x)


def build_table(f, a, b, m_1):
    table = {}
    for i in range(m_1):
        x = a + i * (b - a) / m_1
        table[x] = f(x)
    return table


def lagrange_polynomial(table, x):
    result = 0
    for i in range(len(table)):
        x_i, y_i = table[i]
        li = 1
        for j in range(len(table)):
            if j != i:
                x_j, _ = table[j]
                li *= (x - x_j) / (x_i - x_j)
        result += y_i * li
    return result


def print_table(table, title="Таблица значений"):
    print(f"\n=== {title} ===")
    print(f"{'x':>10} | {'f(x)':>10}")
    print("-" * 25)
    for x, y in table:
        print(f"{x:10.7f} | {y:10.7f}")
    print()


def main():
    print("=== ЗАДАЧА АЛГЕБРАИЧЕСКОГО ИНТЕРПОЛИРОВАНИЯ === ")
    print("Интерполяционный многочлен в форме Лагранжа")
    print("Функция: f(x) = 1 - exp(-2 * x)")
    print("Интервал: [0, 1]")
    print("n = 7, m + 1 = 15")

    m_1 = int(input("Введите число значений в таблице (m+1): "))
    a = float(input("Введите левую границу отрезка [a, b]: "))
    b = float(input("Введите правую границу отрезка [a, b]: "))

    # Построение таблицы
    table_dict = build_table(f, a, b, m_1)
    table = list(table_dict.items())
    print_table(table, "Исходная таблица значений")

    x = float(input("\nВведите точку интерполирования x: "))
    n = int(
        input(f"Введите степень интерполяционного многочлена n ≤ {m_1 - 1}: "))

    modified_table = sorted(table, key=lambda z: abs(x - z[0]))[:n + 1]
    modified_table = sorted(modified_table, key=lambda z: z[0])
    print_table(modified_table, "Обновленная таблица значений")

    y_interp = lagrange_polynomial(modified_table, x)
    print(
        f"\nЗначение интерполяционного многочлена в точке x = {x}: {y_interp}")
    print(f"Точное значение функции в точке x = {x}: {f(x)}")
    print(f"Погрешность интерполяции: {abs(y_interp - f(x))}")


if __name__ == "__main__":
    main()
