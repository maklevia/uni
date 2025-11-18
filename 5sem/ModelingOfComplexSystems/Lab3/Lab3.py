import numpy as np
import csv
from datetime import datetime
import time
import matplotlib.pyplot as plt
import os
import pandas as pd


# -------------------- Підготовка директорії --------------------
def ensure_logs_directory():
    if not os.path.exists('logs'):
        os.makedirs('logs')


# -------------------- Читання даних --------------------
def read_file(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
        input_data = []
        for line in lines:
            values = line.strip().split()
            input_data.append([float(v) for v in values])
    return np.array(input_data).T  # повертаємо у формат (N,6)


# -------------------- Матриця системи --------------------
def init_matr(params):
    c1, c2, c3, c4 = params['c1'], params['c2'], params['c3'], params['c4']
    m1, m2, m3 = params['m1'], params['m2'], params['m3']
    A = np.array([
        [0, 1, 0, 0, 0, 0],
        [-(c1 + c2) / m1, 0, c2 / m1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [c2 / m2, 0, -(c2 + c3) / m2, 0, c3 / m2, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, c3 / m3, 0, -(c3 + c4) / m3, 0]
    ])
    return A


# -------------------- RK4 для y та чутливостей --------------------
def get_y(a_matr, y_cur, h):
    k1 = h * np.dot(a_matr, y_cur)
    k2 = h * np.dot(a_matr, y_cur + k1 / 2)
    k3 = h * np.dot(a_matr, y_cur + k2 / 2)
    k4 = h * np.dot(a_matr, y_cur + k3)
    return y_cur + (k1 + 2 * k2 + 2 * k3 + k4) / 6


def get_u_matr(a_matr, b_matr, u_matr, h):
    k1 = h * (np.dot(a_matr, u_matr) + b_matr)
    k2 = h * (np.dot(a_matr, u_matr + k1 / 2) + b_matr)
    k3 = h * (np.dot(a_matr, u_matr + k2 / 2) + b_matr)
    k4 = h * (np.dot(a_matr, u_matr + k3) + b_matr)
    return u_matr + (k1 + 2 * k2 + 2 * k3 + k4) / 6


# -------------------- Finite differences --------------------
def finite_diff(y_vec_func, beta_keys, beta_values, delta=1e-5):
    n = len(y_vec_func(beta_values))
    m = len(beta_keys)
    deriv_matrix = np.zeros((n, m))

    for j in range(m):
        key = beta_keys[j]
        original_value = beta_values[key]
        beta_values[key] = original_value + delta
        y_plus = y_vec_func(beta_values)
        beta_values[key] = original_value - delta
        y_minus = y_vec_func(beta_values)
        beta_values[key] = original_value
        deriv_matrix[:, j] = (y_plus - y_minus) / (2 * delta)

    return deriv_matrix


# -------------------- Модель рішення --------------------
def get_model_solution(params, y0, t_points, h=0.2):
    a_matrix = init_matr(params)
    y_cur = y0
    y_sol = [y0]
    for _ in range(len(t_points) - 1):
        y_cur = get_y(a_matrix, y_cur, h)
        y_sol.append(y_cur)
    return np.array(y_sol)


# -------------------- Основний алгоритм апроксимації --------------------
def approximate(y_matr, params_fixed, beta_symbols, beta_values, eps=1e-6, h=0.2):
    ensure_logs_directory()
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Для збереження ітерацій
    iterations_list = []

    iteration = 0
    beta_vector = np.array([beta_values[key] for key in beta_symbols])

    def y_vec_func(b_vals_dict):
        all_params = {**params_fixed, **b_vals_dict}
        A = init_matr(all_params)
        return A @ y_approx

    while True:
        all_params = {**params_fixed, **beta_values}
        A_complete = init_matr(all_params)
        y_approx = y_matr[0]
        u_matr = np.zeros((6, len(beta_symbols)))
        integral_inv = np.zeros((len(beta_symbols), len(beta_symbols)))
        integral_mult = np.zeros(len(beta_symbols))
        quality = 0.0

        for i in range(len(y_matr)):
            b_derivative = finite_diff(y_vec_func, beta_symbols, beta_values)
            integral_inv += u_matr.T @ u_matr
            integral_mult += u_matr.T @ (y_matr[i] - y_approx)
            quality += np.sum((y_matr[i] - y_approx) ** 2)
            u_matr = get_u_matr(A_complete, b_derivative, u_matr, h)
            y_approx = get_y(A_complete, y_approx, h)

        integral_inv *= h
        integral_mult *= h
        quality *= h

        delta_beta = np.linalg.inv(integral_inv) @ integral_mult
        beta_vector += delta_beta
        beta_values = {beta_symbols[i]: beta_vector[i] for i in range(len(beta_symbols))}

        # Зберігаємо ітерацію
        iter_record = {'Iteration': iteration}
        for i, key in enumerate(beta_symbols):
            iter_record[key] = beta_vector[i]
        iter_record['Quality'] = quality
        iterations_list.append(iter_record)

        iteration += 1
        if quality < eps:
            break

    execution_time = time.time() - start_time

    # Таблиця ітерацій
    results_df = pd.DataFrame(iterations_list)
    print("\nТаблиця ітерацій:")
    print(results_df.to_string(index=False))

    return beta_values, iteration, execution_time, quality, results_df



def plot_three_signals(t_points, measured_data, model_data, save_path=None):
    """
    Побудова трьох графіків на одному полотні для сигналів y1, y2, y3.

    Args:
        t_points (array): масив часу
        measured_data (array): вхідні дані (N,6) або (N,3) — беремо перші 3 сигнали
        model_data (array): модельні сигнали (N,6) або (N,3)
        save_path (str, optional): якщо вказано — зберігає графік у файл
    """
    labels = ['y1', 'y2', 'y3']
    plt.figure(figsize=(10, 8))

    for i in range(3):
        plt.subplot(3, 1, i + 1)
        plt.plot(t_points, measured_data[:, i], 'r.', label='Вхідні дані', markersize=4)
        plt.plot(t_points, model_data[:, i], 'b-', label='Модель', linewidth=1.5)
        plt.ylabel(labels[i])
        plt.legend()
        plt.grid(True)

    plt.xlabel('Час t')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()


# -------------------- Головна функція --------------------
def main():
    ensure_logs_directory()

    # Данні
    y_data = read_file('./y1.txt')  # (N,6)
    t_points = np.arange(len(y_data)) * 0.2

    # Фіксовані параметри
    params_fixed = {'c2': 0.3, 'c3': 0.2, 'c4': 0.12, 'm3': 18.0}  # c3 фіксоване
    beta_symbols = ['c1', 'm1', 'm2']  # шукаються тільки ці
    beta_values = {'c1': 0.1, 'm1': 11.0, 'm2': 23.0}

    result, iters, exec_time, quality, results_df = approximate(y_data, params_fixed, beta_symbols, beta_values, 1e-6)

    final_params = {**params_fixed, **result}
    model_solution = get_model_solution(final_params, y_data[0], t_points)

    print("\nIdentified parameters:")
    for k, v in result.items():  # виводимо тільки c1, m1, m2
        print(f"{k}: {v:.6f}")
    print(f"Quality indicator: {quality:.6e}")
    # measured_data і model_data — це масиви розміром (N,6), беремо перші 3 сигнали
    plot_three_signals(t_points, y_data, model_solution)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


if __name__ == "__main__":
    main()
