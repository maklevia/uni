import numpy as np

# Перевірка матриці на невиродженість
def is_matrix_nonsingular(matrix):
    det = np.linalg.det(matrix)
    return abs(det) > 0

# Перевірка матриці на симетричність
def is_matrix_symmetric(matrix):
    return np.allclose(matrix, matrix.T)

# Перевірка матриці на додатньо визначеність
def is_positive_definite(matrix):
    n = matrix.shape[0]
    for i in range(1, n + 1):
        minor = matrix[:i, :i]
        det = np.linalg.det(minor)
        if det <= 0:
            return False
    return True

# Обчислення норми вектора (максимальне значення компоненти)
def vector_norm(vector):
    return np.max(np.abs(vector))

# Обчислення норми матриці (максимальна сума модулів елементів по рядках)
def matrix_norm(matrix):
    return np.max(np.sum(np.abs(matrix), axis=1))

# Обчислення числа обумовленості матриці
def condition_number(A):
    try:
        # Обернена матриця
        A_inv = np.linalg.inv(A)
        
        # Норми матриць
        norm_A = matrix_norm(A)
        norm_A_inv = matrix_norm(A_inv)
        
        # Число обумовленості
        cond_number = norm_A * norm_A_inv
        
        return cond_number, A_inv, norm_A, norm_A_inv
        
    except np.linalg.LinAlgError:
        print("Помилка: матриця вироджена, оберненої матриці не існує")
        return None, None, None, None


# Метод квадратних коренів
def square_root_method(A, b, n):
    # Перевірка матриці А на симетричність
    is_symmetric = is_matrix_symmetric(A)
    if is_symmetric:
        print("Матриця А симетрична. Метод квадратних коренів можна застосувати.")
    else:
        print("Матриця А не симетрична. Метод квадратних коренів використати не можна.")
        return None

    # Ініціалізуємо матриці S та D
    S = np.zeros((n, n))
    D = np.zeros((n, n))
    
    # Обчислюємо елементи матриць S та D
    for i in range(n):
        # Обчислюємо суму для діагонального елемента
        sum_diag = 0
        for k in range(i):
            sum_diag += (S[k, i] ** 2) * D[k, k]
        
        # Обчислюємо d_ii
        dii_value = A[i, i] - sum_diag
        D[i, i] = np.sign(dii_value)
        
        # Обчислюємо s_ii
        S[i, i] = np.sqrt(abs(dii_value))
        
        # Обчислюємо недіагональні елементи s_ij
        for j in range(i + 1, n):
            sum_off_diag = 0
            for p in range(i):
                sum_off_diag += S[p, i] * D[p, p] * S[p, j]
            
            if D[i, i] * S[i, i] != 0:
                S[i, j] = (A[i, j] - sum_off_diag) / (D[i, i] * S[i, i])
            else:
                S[i, j] = 0
    
    print("Матриця S:")
    print(np.round(S, 3))
    print("\nМатриця D:")
    print(np.round(D, 3))
    print()
    
    # Обчислюємо визначник матриці
    det_A = 1.0
    for k in range(n):
        det_A *= D[k, k] * (S[k, k] ** 2)
    print(f"Визначник матриці A = {det_A:.6f}")

    
    # Розв'язуємо систему S^T * D * y = b
    ST_D = S.T @ D

    y = np.zeros(n)
    for i in range(n):
        sum_val = 0
        for j in range(i):
            sum_val += ST_D[i, j] * y[j]
        y[i] = (b[i] - sum_val) / ST_D[i, i]
    
    print("Вектор y:")
    print(np.round(y, 3))
    print()
    
    # Розв'язуємо систему S * x = y
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        sum_val = 0
        for j in range(i + 1, n):
            sum_val += S[i, j] * x[j]
        x[i] = (y[i] - sum_val) / S[i, i]
    
    return x, S, D


# Метод Зейделя для розв'язання СЛАР
def seidel_method(A, b, epsilon, n, max_iterations=1000):
    # Перевірка умов збіжності
    is_pos_def = is_positive_definite(A)
    if is_pos_def:
        print("Умова збіжності (за додатньо визначеною матрицею А) виконується\n")
    else:
        print("Умови збіжності (за додатньо визначеною матрицею А) не виконуються!")
        return None, 0, False

    # Початкове наближення (нульовий вектор)
    x_old = np.zeros(n)
    x_new = np.zeros(n)

    print("Ітераційний процес:")
    for iteration in range(max_iterations):
        x_old = x_new.copy()

        for i in range(n):
            sum1 = 0  # Сума для j = 1 до i-1 (використовуємо x_new)
            for j in range(i):
                sum1 += (A[i, j] / A[i, i]) * x_new[j]

            sum2 = 0  # Сума для j = i+1 до n (використовуємо x_old)
            for j in range(i + 1, n):
                sum2 += (A[i, j] / A[i, i]) * x_old[j]

            x_new[i] = -sum1 - sum2 + b[i] / A[i, i]

        # Обчислюємо норму різниці
        diff_norm = vector_norm(x_new - x_old)

        print(f"Ітерація {iteration + 1}: x = {np.round(x_new, 5)}, ||x^(k+1) - x^(k)|| = {diff_norm:.6f}")

        # Перевіряємо умову припинення
        if diff_norm <= epsilon:
            print(f"\nДосягнуто задану точність ε = {epsilon}")
            print(f"Кількість ітерацій: {iteration + 1}")
            return x_new, iteration + 1, True

    print(f"\nДосягнуто максимальну кількість ітерацій: {max_iterations}")
    return x_new, max_iterations, False




def main():
    # Задаємо симетричну матрицю 4x4 з коефіцієнтами
    A = np.array([
        [4, 1, 2, 1],
        [1, 5, 0, 2],
        [2, 0, 22, 1],
        [1, 2, 1, 12]
    ], dtype=float)
    
    # Задаємо вектор правої частини СЛАР
    b = np.array([10, 8, 12, 15], dtype=float)

    # Розмірність матриці А
    n = A.shape[0]
    
    print("Матриця коефіцієнтів A:")
    print(np.round(A, 3))
    print("\nВектор правої частини b:")
    print(np.round(b, 3))
    print()

    # Перевіряємо матрицю на невиродженість
    is_nonsingular = is_matrix_nonsingular(A)
    if is_nonsingular:
        print("Матриця невироджена")
        
        # Обчислюємо число обумовленості
        print("\nОбчислення числа обумовленості:")
        cond_num, A_inv, norm_A, norm_A_inv = condition_number(A)
        
        if cond_num is not None:
            print("Обернена матриця A^(-1):")
            print(np.round(A_inv, 3))
            print(f"\nНорма початкової матриці A: ||A|| = {norm_A:.6f}")
            print(f"Норма оберненої матриці A^(-1): ||A^(-1)|| = {norm_A_inv:.6f}")
            print(f"Число обумовленості: cond(A) = ||A|| * ||A^(-1)|| = {cond_num:.6f}")
            
            # Інтерпретація числа обумовленості
            if cond_num < 1e12:
                print("Матриця добре обумовлена")
            else:
                print("Матриця погано обумовлена")
        print()
        
    else:
        print("Матриця вироджена")
        return

    # Розв'язуємо СЛАР методом квадратного кореня
    print("\nРозв'язання СЛАР методом квадратного кореня:")
    result = square_root_method(A, b, n)
    
    if result is not None:
        x, S, D = result
        print("Розв'язок x:")
        print(np.round(x, 3))
        print()

        # Перевіряємо розв'язок
        print("Перевірка розв'язку (A * x):")
        print(np.round(A @ x, 3))
        print()
    
    # Розв'язуємо СЛАР методом Зейделя
    print("Розв'язання СЛАР методом Зейделя:")
    epsilon = 0.5
    x_seidel, iterations, converged = seidel_method(A, b, epsilon, n)
    
    if converged:
        print(f"\nРозв'язок методом Зейделя (ε = {epsilon}):")
        print(np.round(x_seidel, 5))
        
        # Перевіряємо розв'язок
        print("\nПеревірка розв'язку (A * x):")
        print(np.round(A @ x_seidel, 5))
    else:
        print("Метод Зейделя не збігся за задану кількість ітерацій")




if __name__ == "__main__":
    main()
