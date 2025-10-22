import numpy as np

# Перевірка матриці на симетричність (для методу обертань Якобі)
def is_symmetric(A, tolerance=1e-10):
    return np.allclose(A, A.T, atol=tolerance)

# Метод обертань (Якобі)
def jacobi_method(A, epsilon, max_iterations=1500):
    
    # Перевірка симетричності
    if not is_symmetric(A):
        raise ValueError("Матриця не є симетричною! Метод Якобі можна застосовувати тільки до симетричних матриць.")
    
    # Ініціалізація
    n = A.shape[0]
    A_k = A.copy()
    U = np.eye(n)  # Матриця власних векторів
    
    for iteration in range(max_iterations):
        # Знаходимо максимальний за модулем позадіагональний елемент
        max_val = 0
        i_max, j_max = 0, 0
        
        for i in range(n):
            for j in range(i + 1, n):  # Тільки верхня трикутна частина
                if abs(A_k[i, j]) > max_val:
                    max_val = abs(A_k[i, j])
                    i_max, j_max = i, j
        
        # Перевірка умови припинення
        # Обчислюємо суму квадратів позадіагональних елементів
        off_diagonal_sum = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    off_diagonal_sum += A_k[i, j] ** 2
        
        if off_diagonal_sum <= epsilon:
            print(f"Досягнуто точність після {iteration + 1} ітерацій")
            break
        
        # Обчислення кута обертання
        if abs(A_k[i_max, i_max] - A_k[j_max, j_max]) < 1e-15:
            phi = np.pi / 4  # Якщо діагональні елементи рівні
        else:
            phi = 0.5 * np.arctan(2 * A_k[i_max, j_max] / (A_k[i_max, i_max] - A_k[j_max, j_max]))
        
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        # Створення матриці обертання U_k
        U_k = np.eye(n)
        U_k[i_max, i_max] = cos_phi
        U_k[i_max, j_max] = sin_phi
        U_k[j_max, i_max] = -sin_phi
        U_k[j_max, j_max] = cos_phi
        
        # Обертання матриці: A_{k+1} = U_k^T * A_k * U_k
        A_k = U_k.T @ A_k @ U_k
        
        # Накопичення власних векторів: U = U * U_k
        U = U @ U_k

    # Власні значення - діагональні елементи останньої матриці
    eigenvalues = np.diag(A_k)

    # Сортування за спаданням модуля
    idx = np.argsort(np.abs(eigenvalues))[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = U[:, idx]

    return eigenvalues, eigenvectors, iteration + 1

# Степеневий метод для пошуку найбільшого власного значення
def power_method_max_eigenvalue(A, x0, epsilon, max_iterations=1000):
    x_k = x0.copy()
    lambda_k_prev = 0

    for iteration in range(max_iterations):
        # x_{k+1} = A * x_k
        x_k_plus_1 = A @ x_k

        # Обчислення lambda_k
        lambda_k = x_k_plus_1[0] / x_k[0]

        # Нормування вектора (ділимо на найбільшу компоненту)
        max_component = np.max(np.abs(x_k_plus_1))
        e = x_k_plus_1 / max_component


        print(f"Ітерація {iteration + 1}:")
        print(f"  e_{iteration + 1} = {np.round(e, 6)}")
        print(f"  lambda_{iteration + 1} = {lambda_k:.6f}")

        # Перевірка умови припинення
        if iteration > 0 and abs(lambda_k - lambda_k_prev) <= epsilon:
            print(f"\nДосягнуто точність після {iteration + 1} ітерацій")
            print(f"|lambda_{iteration + 1} - lambda_{iteration}| = {abs(lambda_k - lambda_k_prev):.6f} <= {epsilon}")
            break

        lambda_k_prev = lambda_k
        x_k = e

    return lambda_k, x_k, iteration + 1

# Пошук найменшого власного значення
def power_method_min_eigenvalue(A, x0, epsilon=0.01, max_iterations=1000):

    # Спочатку знаходимо найбільше власне значення
    n = A.shape[0]
    lambda_max, x_max, iterations_max = power_method_max_eigenvalue(A, x0, epsilon)

    # Будуємо матрицю B = lambda_max * E - A
    n = A.shape[0]
    E = np.eye(n)
    B = lambda_max * E - A

    print(f"\nМатриця B = lambda_max * E - A:")
    print(np.round(B, 6))

    # Знаходимо найбільше власне значення матриці B
    print("Ітерації для пошуку найбільшого власного значення матриці B")
    lambda_max_B, x_max_B, iterations_min = power_method_max_eigenvalue(B, x0, epsilon)

    # Найменше власне значення A
    lambda_min = lambda_max - lambda_max_B
    print(f"Найменше власне значення матриці A: {lambda_min:.6f}")

    return lambda_min, lambda_max, iterations_max, iterations_min


def main():
    # Задаємо симетричну матрицю 4x4 з коефіцієнтами
    A = np.array([
        [4, 1, 2, 1],
        [1, 5, 0, 2],
        [2, 0, 12, 1],
        [1, 2, 1, 22]
    ], dtype=float)

    print("Матриця A:")
    print(np.round(A, 3))

    # Параметри методу
    epsilon = 0.01
    x0 = np.array([1, 1, 1, 1], dtype=float)

    print(f"Точність (epsilon): {epsilon}")
    print(f"Початкове наближення: {x0}")

    # Розмірність матриці А
    n = A.shape[0]

    # Знаходимо найбільше власне значення
    print("Степеневий метод для найбільшого власного значення:")
    lambda_max, x_max, iterations_max = power_method_max_eigenvalue(A, x0, epsilon)

    # Знаходимо найменше власне значення
    lambda_min, lambda_max_check, iterations_max_check, iterations_min = power_method_min_eigenvalue(A, x0, epsilon)

    # Підсумок результатів степеневого методу
    print("\nРезультати степеневого методу")
    print(f"Найбільше власне значення: {lambda_max:.6f}")
    print(f"Найменше власне значення: {lambda_min:.6f}")
    print(f"Кількість ітерацій для lambda_max: {iterations_max}")
    print(f"Кількість ітерацій для lambda_min: {iterations_min}")

    print("\nМетод обертань (Якобі)")

    try:
        # Метод Якобі
        eigenvalues_jacobi, eigenvectors_jacobi, iterations_jacobi = jacobi_method(A, epsilon)
        print(f"Власні значення: {np.round(eigenvalues_jacobi, 6)}")
        print(f"Кількість ітерацій: {iterations_jacobi}")

        print(f"\nМатриця власних векторів U:")
        print(np.round(eigenvectors_jacobi, 6))

    except ValueError as e:
        print(f"Помилка: {e}")
        print("Метод Якобі не можна застосувати до цієї матриці.")


    print("\nПорівняння результатів методів")
    if 'eigenvalues_jacobi' in locals():
        print(f"\nПорівняння методів:")
        print(f"Метод Якобі - всі власні значення: {np.round(eigenvalues_jacobi, 6)}")
        print(f"Степеневий метод - найбільше: {lambda_max:.6f}, найменше: {lambda_min:.6f}")



if __name__ == "__main__":
    main()