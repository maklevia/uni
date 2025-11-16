import numpy as np
import matplotlib.pyplot as plt


# Функція
def f(x):
    return x**4 - 5.74*x**3 + 8.18*x - 3.48*np.cos(x)


# --------------------- ЛАГРАНЖ ---------------------
# Класичний поліном Лагранжа з методички
def lagrange_polynomial(x, nodes, values):
    n = len(nodes)
    P = np.zeros_like(x)

    for i in range(n):
        Li = np.ones_like(x)
        for j in range(n):
            if i != j:
                Li *= (x - nodes[j]) / (nodes[i] - nodes[j])
        P += values[i] * Li

    return P


# --------------------- ВУЗЛИ ЧЕБИШЕВА ---------------------
# Формула з методички: x_k = (a+b)/2 + (b−a)/2·cos((2k+1)π/(2n))
def chebyshev_nodes(a, b, n):
    nodes = []
    for k in range(n):
        xk = (a + b)/2 + (b - a)/2 * np.cos((2*k + 1)*np.pi / (2*n))
        nodes.append(xk)
    return np.array(sorted(nodes))


# --------------------- ПРЯМА ІНТЕРПОЛЯЦІЯ ---------------------
# Розв'язання L(x)=0 методом Ньютона (дозволяється в методичці)
def direct_interpolation(nodes, values, eps=1e-10):

    def L(x):
        return lagrange_polynomial(np.array([x]), nodes, values)[0]

    def dL(x):
        h = 1e-7
        return (L(x+h) - L(x-h)) / (2*h)

    x = np.mean(nodes)

    while True:
        x_next = x - L(x) / dL(x)
        if abs(x - x_next) < eps:
            return x_next
        x = x_next


# --------------------- ПЕРЕВІРКА МОНОТОННОСТІ ---------------------
def is_monotonic(a, b, func, n_samples=1000, tol=1e-10):
    """
    Перевіряє, чи є функція монотонною на проміжку [a, b]
    
    Параметри:
    -----------
    a, b : float
        Межі проміжку
    func : callable
        Функція для перевірки
    n_samples : int
        Кількість точок для перевірки
    tol : float
        Допустима похибка для визначення знаку похідної
    
    Повертає:
    ---------
    tuple : (bool, str)
        (is_monotonic, type) де type може бути 'increasing', 'decreasing', або 'non_monotonic'
    """
    x_samples = np.linspace(a, b, n_samples)
    y_samples = func(x_samples)
    
    # Обчислюємо похідну в кожній точці (чисельне диференціювання)
    h = (b - a) / n_samples
    derivatives = np.gradient(y_samples, h)
    
    # Перевіряємо знак похідної
    positive_count = np.sum(derivatives > tol)
    negative_count = np.sum(derivatives < -tol)
    zero_count = np.sum(np.abs(derivatives) <= tol)
    
    # Якщо всі похідні додатні (або нульові) - функція зростає
    if negative_count == 0:
        return True, 'increasing'
    # Якщо всі похідні від'ємні (або нульові) - функція спадає
    elif positive_count == 0:
        return True, 'decreasing'
    # Якщо є і додатні, і від'ємні - функція не монотонна
    else:
        return False, 'non_monotonic'


# --------------------- ОБЕРНЕНА ІНТЕРПОЛЯЦІЯ ---------------------
# Формула з методички:
# x(y) = Σ x_i * ( ω(y)/((y−y_i)·ω'(y_i)) )
def inverse_interpolation(nodes, values, y_target=0):

    n = len(nodes)

    def omega(y):
        w = 1
        for yi in values:
            w *= (y - yi)
        return w

    def omega_prime(i):
        w = 1
        for j in range(n):
            if j != i:
                w *= (values[i] - values[j])
        return w

    w_y = omega(y_target)
    x_res = 0

    for i in range(n):
        x_res += nodes[i] * w_y / ((y_target - values[i]) * omega_prime(i))

    return x_res


# --------------------- ВИВІД ПОЛІНОМА ЛАГРАНЖА ---------------------
def print_lagrange_polynomial(nodes, values):
    """
    Виводить поліном Лагранжа в символьному вигляді
    P(x) = Σ(i=0 to n) f(x_i) * L_i(x)
    """
    n = len(nodes)
    print("\n" + "=" * 80)
    print("ПОЛІНОМ ЛАГРАНЖА")
    print("=" * 80)
    print(f"\nP(x) = Σ(i=0 to {n-1}) f(x_i) * L_i(x)\n")

    terms = []
    for i in range(n):
        # Формуємо L_i(x)
        Li_parts = []
        denominator = 1.0

        for j in range(n):
            if i != j:
                Li_parts.append(f"(x - {nodes[j]:.8f})")
                denominator *= (nodes[i] - nodes[j])

        Li_str = " * ".join(Li_parts) if Li_parts else "1"

        # Формуємо повний член
        coeff = values[i] / denominator
        term_str = f"f(x_{i}) * L_{i}(x)"
        terms.append((i, term_str, coeff, Li_str, denominator))

        print(f"L_{i}(x) = ", end="")
        if len(Li_parts) > 0:
            print(f"({' * '.join(Li_parts)}) / {denominator:.8e}")
        else:
            print("1")

    print("\n" + "-" * 80)
    print("Повний поліном:")
    print("-" * 80)
    print("\nP(x) = ", end="")

    for i, (idx, term_str, coeff, Li_str, denom) in enumerate(terms):
        if i > 0:
            print(" + ", end="")
        print(f"({values[idx]:.8e}) * L_{idx}(x)", end="")

    print("\n\nРозкритий вигляд:")
    print("-" * 80)
    for i, (idx, term_str, coeff, Li_str, denom) in enumerate(terms):
        if i > 0:
            print(" + ", end="")
        else:
            print("P(x) = ", end="")

        if len(Li_str) > 1 and Li_str != "1":
            print(f"({values[idx]:.8e}) * ({Li_str}) / ({denom:.8e})", end="")
        else:
            print(f"{values[idx]:.8e}", end="")

    print("\n" + "=" * 80 + "\n")


# --------------------- ГРАФІК ---------------------
def draw(a, b, nodes=None, values=None, ylim=None):
    x = np.linspace(a, b, 2000)
    y = f(x)

    plt.figure(figsize=(12, 7))
    plt.plot(x, y, label="f(x)")

    if nodes is not None:
        xL = np.linspace(nodes[0], nodes[-1], 2000)
        yL = lagrange_polynomial(xL, nodes, values)
        plt.plot(xL, yL, '--', label="Поліном Лагранжа")
        plt.scatter(nodes, values, c='red')

    plt.xlim(a, b)
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])

    plt.grid()
    plt.legend()
    plt.show()


# --------------------- MAIN ---------------------
def main():

    a, b = -2, -1   # згідно методички — інтервал для дослідження кореня
    n = 10             # кількість вузлів Чебишева

    print("Вузли Чебишева:")
    nodes = chebyshev_nodes(a, b, n)
    values = f(nodes)
    for i in range(n):
        print(f"x[{i}] = {nodes[i]:.8f},  f(x) = {values[i]:.8f}")

    # Виведення полінома Лагранжа
    print_lagrange_polynomial(nodes, values)

    draw(-10, 10, nodes, values, ylim=(-5, 50))
    draw(-1.2, -1, nodes, values, ylim=(0, -2.5))

    print("\nПряма інтерполяція...")
    x_direct = direct_interpolation(nodes, values)
    print(f"x* = {x_direct:.10f},   f(x*) = {f(x_direct):.2e}")

    # Перевірка монотонності перед оберненою інтерполяцією
    print("\nПеревірка монотонності функції на проміжку...")
    is_mono, mono_type = is_monotonic(a, b, f)
    
    if mono_type == 'increasing':
        print(f"Функція монотонно зростає на проміжку [{a}, {b}]")
    elif mono_type == 'decreasing':
        print(f"Функція монотонно спадає на проміжку [{a}, {b}]")
    else:
        print(f"Функція НЕ є монотонною на проміжку [{a}, {b}]")
    
    print("\nОбернена інтерполяція...")
    if is_mono:
        x_inverse = inverse_interpolation(nodes, values, 0)
        print(f"x* = {x_inverse:.10f},   f(x*) = {f(x_inverse):.2e}")
        
        print("\nПорівняння:")
        print(f"Різниця |Δ| = {abs(x_direct - x_inverse):.3e}")
    else:
        print("️  Обернена інтерполяція НЕ використовується, оскільки функція не є монотонною на проміжку!")


if __name__ == "__main__":
    main()
