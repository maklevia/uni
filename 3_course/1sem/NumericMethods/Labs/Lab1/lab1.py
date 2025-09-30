import numpy as np
import matplotlib.pyplot as plt


# Означення функції та її похідних
def f(x):
    return x ** 3 - 3 * x ** 2 - 17 * x + 22


def df(x):
    return 3 * x ** 2 - 6 * x - 17


def d2f(x):
    return 6 * x - 6


# Побудова графіка функції та її похідних
def draw(a, b):
    x = np.linspace(0, 3, 1000)
    y1 = f(x)
    y2 = df(x)
    y3 = d2f(x)

    fig, ax = plt.subplots()
    ax.plot(x, y1, label='f(x)=x^3-3x^2-17x+22')
    ax.plot(x, y2, label="f'(x)=3x^2-6x-17")
    ax.plot(x, y3, label="f''(x)=6x-6")

    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    # Вертикальні лінії x=a та x=b (пунктир, червоний)
    ax.axvline(a, color='red', linestyle='--', linewidth=1)
    ax.axvline(b, color='red', linestyle='--', linewidth=1)

    ax.set_title('Plotting Functions', size=14)
    ax.set_xlim(0, 3)
    ax.set_ylim(-20, 20)
    plt.title('Графік функції та її похідних')
    plt.legend()
    plt.show()


# Перевірка наявності кореня на проміжку
def has_root(f, a, b) -> bool:
    if a >= b:
        raise ValueError("The interval must satisfy 'a < b'.")

    if f(a) * f(b) <= 0:
        return True
    else:
        return False


# Знаходження мінімуму та максимуму функції на проміжку
def find_min_max(f, a, b, num_points=10000):
    x_values = np.linspace(a, b, num_points)
    derivatives = f(x_values)
    abs_derivatives = np.abs(derivatives)
    Min = np.min(abs_derivatives)
    Max = np.max(abs_derivatives)
    return Min, Max


# Перевірка неперервності функції на проміжку
def continuity(f, a, b, epsilon):
    if not callable(f):
        raise ValueError("The provided 'f' is not a callable function.")

    def limit(f, c, delta_x=1e-6):
        left_limit = f(c - delta_x)
        right_limit = f(c + delta_x)
        return left_limit, right_limit

    if a <= b and (a is None or b is None):
        return False

    left_limit_a, right_limit_a = limit(f, a)
    left_limit_b, right_limit_b = limit(f, b)

    # Перевірка чи збігаються значення лімітів зі значеннями функції на кінцях відрізка
    if abs(left_limit_a - f(a)) > epsilon or abs(right_limit_a - f(a)) > epsilon or abs(
            left_limit_b - f(b)) > epsilon or abs(right_limit_b - f(b)) > epsilon:
        return False

    return True


# Перевірка знакосталості функції
def check_sign_changes(f, a, b, num_points=10000) -> bool:
    x_values = np.linspace(a, b, num_points)
    signs = [np.sign(f(x)) for x in x_values]

    sign_changes = False
    for i in range(1, len(signs)):
        if signs[i] != signs[i - 1]:
            sign_changes = True
            break

    return sign_changes


# Метод релаксації
def relaxation_method(x0, epsilon, a, b, f):
    if not continuity(f, a, b, epsilon):
        return None, 0
    print("\nМетод релаксації")

    # Знаходження m1, M1, tau (оптимальне), q
    m1, M1 = find_min_max(df, a, b)
    tau = 2 / (M1 + m1)
    q = (M1 - m1) / (M1 + m1)
    print(f"m1 = {m1:.5f}, M1 = {M1:.5f}, τ = {tau:.5f}, q = {q:.5f}")

    # Перевірка умови збіжності
    if m1 > 0 and q < 1:
        print("Достатня умова збіжності виконується.")
    else:
        print("Достатня умова збіжності не виконується.")
        return None, 0

    # Апріорна оцінка
    delta = max(abs(x0 - a), abs(x0 - b))
    n_aprior = int(np.floor(np.log(delta / epsilon) / np.log(1 / q))) + 1
    print("Апріорна оцінка числа ітерацій n =", n_aprior)

    # Ітераційний процес
    x_n = x0
    iteration = 0
    max_iter = 100  # Обмеження максимальної кількості операцій

    print("\nІтерації:")
    print("n\t x_n\t\t f(x_n)\t")
    f_val = f(x_n)  # Виведення 0 ітерації
    print(f"{iteration}\t {x_n:.6f}\t {f_val:.6f}\t")

    while iteration < max_iter:
        iteration += 1
        sign = 1 if df(x_n) < 0 else -1  # Визначення знаку у формулі ітераційного процесу

        x_next = x_n + sign * tau * f(x_n)
        print(f"{iteration}\t {x_next:.6f}\t {f(x_next):.6f}\t")

        # Перевірка умови завершення
        if abs(x_next - x_n) < epsilon:
            x_n = x_next
            break

        x_n = x_next

    # Апостеріорна оцінка n
    n_aposteriori = iteration
    print(f"\nКорінь: x ≈ {x_n:.6f} після {iteration} ітерацій")
    print(f"Апостеріорна оцінка числа ітерацій  = {n_aposteriori}")

    return x_n, iteration


# Метод Ньютона
def newton_method(x0, epsilon, a, b, f, df, d2f):
    print("\nМетод Ньютона")

    # Перевірки на збіжність
    if check_sign_changes(d2f, a, b):
        print("Друга похідна не є знакосталою на проміжку [a, b]")
        return None, 0

    if f(x0) * d2f(x0) <= 0:
        print("Умова f(x0) * d2f(x0) > 0 не виконана")
        return None, 0

    # Знаходження M2 та m1
    m1, M1 = find_min_max(df, a, b)
    m2, M2 = find_min_max(d2f, a, b)

    delta = max(abs(x0 - a), abs(x0 - b))
    q = (M2 * delta) / (2 * m1)

    print(f"m1 = {m1:.5f}, M2 = {M2:.5f}, q = {q:.5f}")

    if q >= 1:
        print("Умова збіжності q < 1 не виконана")
        return None, 0
    print("Умови збіжності виконано")

    # Апріорна оцінка
    inner_value = (np.log(delta / epsilon) / np.log(1 / q)) + 1
    n_aprior = int(np.floor(np.log2(inner_value))) + 1
    print("Апріорна оцінка числа ітерацій n =", n_aprior)

    # Ітераційний процес
    x_n = x0
    iteration = 0
    max_iter = 100

    print("\nІтерації:")
    print("n\t x_n\t\t f(x_n)\t")
    f_val = f(x_n)
    print(f"{iteration}\t {x_n:.6f}\t {f_val:.6f}\t")

    while iteration < max_iter:
        iteration += 1

        x_next = x_n - f(x_n) / df(x_n)
        print(f"{iteration}\t {x_next:.6f}\t {f(x_next):.6f}\t")

        # Перевірка умови завершення
        if abs(x_next - x_n) < epsilon:
            x_n = x_next
            break

        x_n = x_next

    print(f"\nКорінь: x ≈ {x_n:.6f} після {iteration} ітерацій")

    return x_n, iteration


def main():
    a, b = 1.1, 2
    x0 = 1.1
    epsilon = 0.5

    print("Функція: f(x) = x³ - 3x² - 17x + 22")
    print(f"Проміжок: [a, b] = [{a}, {b}]")
    print(f"Початкова точка: x₀ = {x0}")
    print(f"Точність: ε = {epsilon}")

    draw(a, b)

    # Перевірка належності x0 проміжку [a, b]
    if not (a <= x0 <= b):
        print(f"\nПомилка: x₀ = {x0} не належить проміжку [{a}, {b}].")
        return

    # Перевірка наявності кореня на проміжку [a, b]
    if  not has_root(f, a, b):
        print(f"\nКорінь не існує на проміжку [{a}, {b}]")
        return
    print(f"\nКорінь існує на проміжку [{a}, {b}]")

    # Перевірка неперервності f на [a, b]
    is_continuous = continuity(f, a, b, epsilon)
    if not is_continuous:
        print("Функція не є неперервною на заданому проміжку.")
        return
    print(f"\nНеперервність f на проміжку [{a}, {b}] виконується")

    # Виклик методу релаксації
    root_relax, iters_relax = relaxation_method(x0, epsilon, a, b, f)

    # Виклик методу Ньютона
    root_newton, iters_newton = newton_method(x0, epsilon, a, b, f, df, d2f)

if __name__ == "__main__":
    main()