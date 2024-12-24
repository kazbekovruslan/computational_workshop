import math


def exact_solution(x):
    return (math.sin(x) - math.cos(x) + 1 / math.exp(x)) / 2


def f(x, y):
    return math.sin(x) - y


def derivative(i):
    arr = [0, 0, 1, -1, 0, 0, 1, -1]
    return arr[i]


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
            cell_str = str(cell) if cell is not None else ""
            row_str += f"{cell_str:<{col_widths[i]}} | "
        print(row_str)

    # Печатаем таблицу
    print_separator()  # Верхний разделитель
    print_row(table[0])  # Заголовок таблицы
    print_separator()  # Разделитель под заголовком
    for row in table[1:]:
        print_row(row)  # Строки с данными
    print_separator()  # Нижний разделитель


def taylor_method(x0, y0, h, N):
    x_vals = [x0 + k * h for k in range(-2, N + 1)]
    y_vals = [0] * len(x_vals)
    y_vals[2] = y0
    for k in range(-2, N+1):
        x_k = x0 + k * h
        y_vals[k+2] = sum(derivative(i) * (x_k - x0)**i /
                          math.factorial(i) for i in range(8))
    return x_vals, y_vals


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
    x, y = taylor_method(x0, y0, h, N)

    for i in range(4, len(y)-1):
        y[i+1] = y[i] + h/720 * (1901 * f(x[i], y[i]) -
                                 2774 * f(x[i-1], y[i-1]) +
                                 2616 * f(x[i-2], y[i-2]) -
                                 1274 * f(x[i-3], y[i-3]) +
                                 251 * f(x[i-4], y[i-4]))
    return x, y


def adjust_table(arr, max_len):
    while (len(arr) < max_len):
        arr.insert(0, None)


def main():
    print("=== Задача численного решение задачи Коши для обыкновенного дифференциального уравнения первого порядка ===")
    print("Вариант №5")
    print("ОДУ: y'(x) = -y(x) + sin(x)")
    print("Решение: y(x) = [sin(x) - cos(x) + 1/exp(x)] / 2\n")
    while True:
        try:
            print("Введите параметры задачи:")
            h = float(input("Введите шаг h: "))
            N = int(input("Введите количество шагов N: "))
            x0, y0 = 0, 0

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
            headers = [[
                "Шаг", "x", "Точное значение",
                "Тейлор", "Погрешность Тейлор", "Рунге-Кутта", "Эйлер", "Эйлер I", "Эйлер II", "Адамс 4-го порядка"
            ]]
            table_data += headers
            for i in range(max_len):
                exact = exact_solution(x_t[i])
                table_data.append([
                    i-2, f'{x_t[i]:.4f}', exact,
                    y_taylor[i], str(
                        f'{abs(exact-y_taylor[i]):.13e}'), y_rk4[i], y_euler[i],
                    y_euler1[i], y_euler2[i], y_adams[i]
                ])

            last_idx = max_len - 1
            exact = exact_solution(x_t[last_idx])
            table_data.append([
                "Погрешность",
                None, None,
                str(f'{abs(y_taylor[last_idx] - exact):.13e}'), None,
                str(f'{abs(y_rk4[last_idx] - exact):.13e}'),
                str(f'{abs(y_euler[last_idx] - exact):.13e}'),
                str(f'{abs(y_euler1[last_idx] - exact):.13e}'),
                str(f'{abs(y_euler2[last_idx] - exact):.13e}'),
                str(f'{abs(y_adams[last_idx] - exact):.13e}')
            ])

            print("\nТаблица решений:")
            print_table(table_data)

            choice = input(
                "Введите 1, чтобы выбрать новые параметры N и h; или 2 -- чтобы завершить программу: ")
            if choice == "1":
                print()
                continue
            elif choice == "2":
                print("=== Завершение работы ===")
                break
            else:
                print("Неверный выбор, попробуйте снова.")
        except ValueError as e:
            print(f"Ошибка ввода: {e}. Попробуйте снова.")


if __name__ == "__main__":
    main()
