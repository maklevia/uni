import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def read_data_file(file_path):
    return np.loadtxt(file_path)

# Дискретне перетворення Фур'є
def DFT(signal, N):
    dft = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            dft[k] += signal[n] * np.exp(-2j * np.pi * k * n / N)
    return dft

# Знаходження локального максимуму (суттєва частота)
def find_significant_frequencies(dft, dt, T, N):
    df = 1 / T

    # беремо тільки ліву половину спектра
    f_magnitude = np.abs(dft[:N//2])

    # шукаємо локальні максимуми
    peaks, _ = find_peaks(f_magnitude)

    if len(peaks) > 0:
        k_star = peaks[np.argmax(f_magnitude[peaks])]
        f_star = k_star * df
        return float(f_star)
    else:
        return None


# Знаходження невідомих параметрів
def fit_model(t, observations, peak_frequencies):
    # Лінійна модель: y = a5 + a1 t^3 + a2 t^2 + a3 t + A1 sin(2π f1 t)
    ones_col = np.ones_like(t)
    t3_col = t**3
    t2_col = t**2
    t1_col = t

    if len(peak_frequencies) >= 1 and peak_frequencies[0] is not None:
        f1 = float(peak_frequencies[0])
        sin_col = np.sin(2 * np.pi * f1 * t)
        X = np.column_stack([t3_col, t2_col, t1_col, sin_col, ones_col])
    else:
        X = np.column_stack([t3_col, t2_col, t1_col, ones_col])

    params, *_ = np.linalg.lstsq(X, observations, rcond=None)
    params = np.round(params).astype(int)
    model_values = X @ params
    return params, model_values


# Побудова графіка спектру
def plot_spectrum(spectrum, dt, dominant_freq):
    magnitude = np.abs(spectrum)
    freqs = np.fft.fftfreq(len(spectrum), dt)
    plt.figure()
    half = len(spectrum)//2
    plt.plot(freqs[:half], magnitude[:half])
    if dominant_freq is not None:
        # find nearest index in positive frequencies
        positive_mask = freqs[:half] >= 0
        pos_freqs = freqs[:half][positive_mask]
        pos_magnitude = magnitude[:half][positive_mask]
        idx = (np.abs(pos_freqs - dominant_freq)).argmin()
        plt.scatter([pos_freqs[idx]], [pos_magnitude[idx]], color='red')
    plt.xlabel("Частота [Hz]")
    plt.ylabel("Амплітуда")
    plt.title("Амплітудний спектр")
    plt.grid(True)
    plt.show()

# Побудова графіка порівняння
def plot_comparison(t, signal, model):
    plt.figure()
    plt.plot(t, signal, label="Спостереження")
    plt.plot(t, model, "--", label="Модель")
    plt.xlabel("Час")
    plt.ylabel("y(t)")
    plt.legend()
    plt.title("Порівняння сигналу та моделі")
    plt.grid(True)
    plt.show()



def main():
    # Читаємо значення y(t) з файлу
    observations = read_data_file('f1.txt')

    # Ініціалізація значень
    T = 5
    N = len(observations)
    dt = 0.01
    t = np.arange(0, N * dt, dt)

    # ДПФ, пошук суттєвої частоти
    spectrum = DFT(observations, N)
    dominant_freq = find_significant_frequencies(spectrum, dt, T, N)
    print("Суттєва частота:", dominant_freq)

    # Графік спектра
    plot_spectrum(spectrum, dt, dominant_freq)

    # Параметри, знайдені методом найменших квадратів
    params, model = fit_model(t, observations, [dominant_freq] if dominant_freq is not None else [])
    print("Оцінені параметри:", params)

    # Графік порівняння
    plot_comparison(t, observations, model)


if __name__ == "__main__":
    main()


