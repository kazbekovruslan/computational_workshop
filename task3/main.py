import math

NUMBER_OF_VARIANT = 5


class Function:
    def __init__(self, f, f_diff, f_2diff):
        self.f = f
        self.f_diff = f_diff
        self.f_2diff = f_2diff


def f_2(x):
    return 1 - math.exp(-2 * x)


def df_2(x):
    return 2 * math.exp(-2 * x)


def ddf_2(x):
    return -4 * math.exp(-2 * x)


def f_1(x):
    return math.exp(1.5 * ((NUMBER_OF_VARIANT % 5) + 1) * x)


def df_1(x):
    return (1.5 * ((NUMBER_OF_VARIANT % 5) + 1)) * math.exp(1.5 * ((NUMBER_OF_VARIANT % 5) + 1) * x)


def ddf_1(x):
    return ((1.5 * ((NUMBER_OF_VARIANT % 5) + 1)) ** 2) * math.exp(1.5 * ((NUMBER_OF_VARIANT % 5) + 1) * x)


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


def create_table(a, h, m_1, func):
    table = []
    for i in range(m_1):
        x_i = a + i * h
        f_x_i = func.f(x_i)
        table.append((x_i, f_x_i))
    return table


def first_derivative(table, h, i):  # формулы 6 7 8
    if i == 0:
        return (-3 * table[i][1] + 4 * table[i + 1][1] - table[i + 2][1]) / (2 * h)
    if i == len(table) - 1:
        return (3 * table[i][1] - 4 * table[i - 1][1] + table[i - 2][1]) / (2 * h)
    return (table[i + 1][1] - table[i - 1][1]) / (2 * h)


def new_first_derivative(table, h, i):
    if i == 0:
        return (-25 * table[i][1] + 48 * table[i+1][1] - 36 * table[i+2][1] + 16 * table[i+3][1] - 3 * table[i+4][1])/(12*h)
    if i == 1:
        return (-3 * table[i-1][1] - 10 * table[i][1] + 18 * table[i+1][1] - 6 * table[i+2][1] + table[i+3][1])/(12*h)
    if i == len(table) - 2:
        return (3 * table[i+1][1] + 10 * table[i][1] - 18 * table[i-1][1] + 6 * table[i-2][1] - table[i-3][1])/(12*h)
    if i == len(table) - 1:
        return (25 * table[i][1] - 48 * table[i-1][1] + 36 * table[i-2][1] - 16 * table[i-3][1] + 3 * table[i-4][1])/(12*h)
    return (table[i-2][1] - 8 * table[i-1][1] + 8 * table[i+1][1] - table[i+2][1])/(12 * h)


def second_derivative(table, h, i):
    if i == 0:
        return (2 * table[i][1] - 5 * table[i+1][1] + 4 * table[i+2][1] - table[i+3][1]) / (h ** 2)
    if i == len(table) - 1:
        return (2 * table[i][1] - 5 * table[i-1][1] + 4 * table[i-2][1] - table[i-3][1]) / (h ** 2)
    return (table[i + 1][1] - 2 * table[i][1] + table[i - 1][1]) / (h ** 2)


def calculate(func, m_1, h, a):
    table = create_table(a, h, m_1, func)

    values_table = [["x_i", "f(x_i)"]]
    for i, (x_i, f_x_i) in enumerate(table):
        values_table.append([f"{x_i}", f"{f_x_i}"])

    print("\nТаблица значений функции:")
    print_table(values_table)

    derivatives_table = [
        ["x_i", "f(x_i)", "df O(h^2)", "|df_true - df|", "df O(h^4)",
         "|df_true - df|", "ddf O(h^2)", "|ddf_true - ddf|"]
    ]
    for i, (x_i, f_x_i) in enumerate(table):
        f_diff = first_derivative(table, h, i)
        f_2diff = second_derivative(table, h, i)
        f_n_diff = new_first_derivative(table, h, i)
        derivatives_table.append([
            f"{x_i}",
            f"{f_x_i:.6f}",
            f"{f_diff}",
            f"{abs(func.f_diff(x_i) - f_diff)}",
            f"{f_n_diff}",
            f"{abs(func.f_diff(x_i) - f_n_diff)}",
            f"{f_2diff}",
            f"{abs(func.f_2diff(x_i) - f_2diff):13f}"
        ])

    print("\nТаблица производных:")
    print_table(derivatives_table)


def calculate_runge(func, m_1, h, a):
    table_h = create_table(a, h, m_1, func)
    table_h2 = create_table(a, h/2, m_1 * 2, func)

    values_table = [["x_i", "f(x_i)"]]

    for i, (x_i, f_x_i) in enumerate(table_h):
        values_table.append([f"{x_i}", f"{f_x_i}"])

    print("\nТаблица значений функции:")
    print_table(values_table)

    derivatives_table_df = [
        ["x_i", "f(x_i)", "J(h)", "|df_true - J(h)|", "J(h/2)", "|df_true - J(h/2)|",
         "J", "|df_true - J|"]
    ]

    derivatives_table_ddf = [
        ["x_i", "f(x_i)", "J(h)", "|ddf_true - J(h)|", "J(h/2)", "|ddf_true - J(h/2)|",
         "J", "|ddf_true - J|"]
    ]

    for i, (x_i, f_x_i) in enumerate(table_h):
        j_h = first_derivative(table_h, h, i)
        j_h_2 = first_derivative(table_h2, h/2, i*2)
        j = (4*j_h_2-j_h)/3

        derivatives_table_df.append([
            f"{x_i}",
            f"{f_x_i:.6f}",
            f"{j_h}",
            f"{abs(func.f_diff(x_i) - j_h)}",
            f"{j_h_2}",
            f"{abs(func.f_diff(x_i) - j_h_2)}",
            f"{j}",
            f"{abs(func.f_diff(x_i) - j)}",
        ])

    print("\nТаблица df:")
    print_table(derivatives_table_df)

    for i, (x_i, f_x_i) in enumerate(table_h):
        j_h_d2 = second_derivative(table_h, h, i)
        j_h_2_d2 = second_derivative(table_h2, h/2, i*2)
        j_d2 = (4*j_h_2_d2-j_h_d2)/3

        derivatives_table_ddf.append([
            f"{x_i}",
            f"{f_x_i:.6f}",
            f"{j_h_d2}",
            f"{abs(func.f_2diff(x_i) - j_h_d2)}",
            f"{j_h_2_d2}",
            f"{abs(func.f_2diff(x_i) - j_h_2_d2)}",
            f"{j_d2}",
            f"{abs(func.f_2diff(x_i) - j_d2)}",
        ])

    print("\nТаблица ddf:")
    print_table(derivatives_table_ddf)


def main(inpt=""):
    print("=== Нахождение производных таблично-заданной функции по формулам численного дифференцирования === ")
    print("Вариант №5")
    while True:
        func_number = 0
        while True:
            print("Введите номер функции, для которой будет решаться задача: ")
            print("1: f(x) = exp(1,5 * x)")
            print("2: f(x) = 1 - exp(-2 * x)")
            func_number = input()
            if func_number == "1":
                func = Function(f_1, df_1, ddf_1)
                break
            elif func_number == "2":
                func = Function(f_2, df_2, ddf_2)
                break
            else:
                print("Некорректный ввод! Введите 1 или 2\n")

        while True:
            runge = False

            while True:
                m_1 = int(
                    input("Введите число значений в таблице (в наших обозначениях это m + 1) >= 5: "))
                if m_1 >= 5:
                    break
                print("Некорректный ввод! m + 1 должно быть >= 5\n")

            a = float(input("Введите левую границу a: "))

            while True:
                h = float(input("Введите размер шага h > 0: "))
                if h > 0:
                    break
                print("Некорректный ввод! h должно быть > 0\n")

            while True:
                if not runge:
                    calculate(func, m_1, h, a)
                else:
                    calculate_runge(func, m_1, h, a)

                while True:
                    inpt = input(
                        "\nВыберите действие:\n1. Начать сначала для другой функции и других значений\n2. Рассмотреть уточнение по Рунге\n3. Выйти\n")
                    if inpt in ["1", "2", "3"]:
                        break
                    print("Некорректный ввод! Введите 1, 2 или 3\n")
                if inpt == "1":
                    break
                if inpt == "2":
                    runge = True
                if inpt == "3":
                    break
            if inpt == "1":
                break
            if inpt == "3":
                break
        if inpt == "3":
            print("=== Выполнение завершено ===")
            break


if __name__ == '__main__':
    main()
