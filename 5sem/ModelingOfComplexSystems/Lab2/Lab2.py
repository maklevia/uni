import numpy as np
import matplotlib.pyplot as plt
import imageio.v3 as iio

# Функція для виведення зображення
def show_image(matrix, title):
    plt.figure()
    plt.imshow(matrix.astype(np.uint8), cmap='gray')
    plt.title(title)
    plt.show(block=True)

# Функція для перевірки правильності знаходження псевдооберненої матрицi
def check_pinv_properties(A, A_plus):
    m, n = A.shape
    # Перевірка розмірів добутків
    AA_plus = A @ A_plus
    A_plusA = A_plus @ A
    if AA_plus.shape != (m, m): return False
    if A_plusA.shape != (n, n): return False
    # Умови Мура–Пенроуза
    if not np.allclose(AA_plus @ A, A): return False
    if not np.allclose(A_plusA @ A_plus, A_plus): return False
    if not np.allclose(AA_plus, AA_plus.T): return False
    if not np.allclose(A_plusA, A_plusA.T): return False
    return True

# Алгоритм Гревіля
def greville(A):
    m, n = A.shape
    # Крок 1: Обробка першого рядка a_1^T
    a1 = A[0, :].reshape(-1, 1)
    denom = a1.T @ a1
    if denom[0, 0] != 0:
        A_plus_i = a1 / denom[0, 0]
    else:
        A_plus_i = np.zeros((n, 1))
    A_i = A[0:1, :]

    # Цикл з додаванням рядків
    for i in range(1, m):
        ai_T = A[i:i + 1, :]
        Z_A = np.eye(n) - A_plus_i @ A_i
        ai = ai_T.T
        # Умова для формули Гревіля
        condition = ai.T @ Z_A @ ai
        # Оновлення A_plus
        if condition[0, 0] > 1e-10:
            v = (Z_A @ ai) / condition[0, 0]

            term = (Z_A @ ai @ ai.T @ A_plus_i) / condition[0, 0]
            A_plus_next = np.hstack([A_plus_i - term, v])
        else:
            R_A = A_plus_i @ A_plus_i.T
            den = 1 + ai.T @ R_A @ ai
            v = (R_A @ ai) / den[0, 0]

            term = (R_A @ ai @ ai.T @ A_plus_i) / den[0, 0]
            A_plus_next = np.hstack([A_plus_i - term, v])
        # Оновлення для наступної ітерації
        A_plus_i = A_plus_next
        A_i = np.vstack([A_i, ai_T])

    return A_plus_i

# Алгоритм за означенням Мура-Пенроуза
def moore_penrose(A, eps=1e-8):
    m, _ = A.shape
    delta = 10.0
    A_prev = A.T @ np.linalg.inv(A @ A.T + (delta**2) * np.eye(m))

    while True:
        delta /= 2
        A_cur = A.T @ np.linalg.inv(A @ A.T + (delta**2) * np.eye(m))
        if np.linalg.norm(A_cur - A_prev) < eps:
            return A_cur
        A_prev = A_cur



# Зчитування вхідного та вихідного зображення
X = np.array(iio.imread("x1.bmp"), dtype=np.float64)
Y = np.array(iio.imread("y1.bmp"), dtype=np.float64)

# Відображення вхідних і вихідних даних
show_image(X, "Вхідне зображення X")
show_image(Y, "Вихідне зображення Y")

# Додаємо рядок одиниць (зсув)
X = np.append(X, np.ones((1, X.shape[1])), axis=0)

m, n = X.shape
p = Y.shape[0]

# Виконуємо алгоритм Гревіля
X_grev_plus = greville(X)
if check_pinv_properties(X, X_grev_plus):
    print("Гревіль: псевдообернена матриця коректна.")

V = np.zeros((p, m))
Z = np.eye(m) - X @ X_grev_plus
A_grev = Y @ X_grev_plus + V @ Z.T

Y_grev = np.clip(A_grev @ X, 0, 255)
show_image(Y_grev, "Результат за Гревілем")

# Виконуємо алгоритм за означенням Мура-Пенроуза
X_mp_plus = moore_penrose(X)
if check_pinv_properties(X, X_mp_plus):
    print("Мур-Пенроуз: псевдообернена матриця коректна.")

Z = np.eye(m) - X @ X_mp_plus
A_mp = Y @ X_mp_plus + V @ Z.T

Y_mp = np.clip(A_mp @ X, 0, 255)
show_image(Y_mp, "Результат за Муром–Пенроузом")


# Порівняння результатів
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.imshow(Y.astype(np.uint8), cmap='gray')
plt.title("Оригінальне Y")

plt.subplot(1, 3, 2)
plt.imshow(Y_grev.astype(np.uint8), cmap='gray')
plt.title("Гревіль")

plt.subplot(1, 3, 3)
plt.imshow(Y_mp.astype(np.uint8), cmap='gray')
plt.title("Мур-Пероуз")

plt.show()
