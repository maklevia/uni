// mk_lab_os_1.cpp : Міжпроцесна комунікація через exit code з використанням polling
//
// Завдання: Реалізувати комунікацію між батьківським та дочірнім процесом
// через код завершення дочірнього процесу з використанням опитування (polling)

#include <iostream>
#include <windows.h>
#include <string>
#include <chrono>
#include <thread>

// Функція для дочірнього процесу
// Виконує просте обчислення та повертає exit code
int childProcess()
{
    std::cout << "[Дочірній процес] Початок виконання..." << std::endl;
    
    // Виконуємо просте обчислення (наприклад, сума чисел від 1 до 10)
    int sum = 0;
    for (int i = 1; i <= 10; ++i)
    {
        sum += i;
        std::cout << "[Дочірній процес] Обчислення: сума = " << sum << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(200)); // Симуляція роботи
    }
    
    std::cout << "[Дочірній процес] Обчислення завершено. Результат: " << sum << std::endl;
    
    // Повертаємо exit code (наприклад, 42 або результат обчислення)
    int exitCode = 42;
    std::cout << "[Дочірній процес] Завершення з кодом: " << exitCode << std::endl;
    
    return exitCode;
}

// Головна функція програми
// Перевіряє аргументи командного рядка для визначення, чи це дочірній процес
int main(int argc, char* argv[])
{
    // Якщо передано аргумент "child", виконуємо код дочірнього процесу
    if (argc > 1 && std::string(argv[1]) == "child")
    {
        int exitCode = childProcess();
        exit(exitCode); // Завершуємо дочірній процес з кодом exitCode
    }
    
    // Інакше виконуємо код батьківського процесу
    std::cout << "=== Міжпроцесна комунікація через Exit Code (Polling) ===" << std::endl;
    std::cout << "[Батьківський процес] Початок роботи..." << std::endl;
    
    // Крок 1: Створення дочірнього процесу
    std::cout << "\n[Батьківський процес] Крок 1: Створення дочірнього процесу..." << std::endl;
    
    // Отримуємо шлях до поточного виконуваного файлу
    char currentExePath[MAX_PATH];
    GetModuleFileNameA(NULL, currentExePath, MAX_PATH);
    
    // Створюємо командний рядок для дочірнього процесу
    // Передаємо аргумент, який вказує, що це дочірній процес
    std::string commandLine = std::string(currentExePath) + " child";
    
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    
    // Ініціалізація структур
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    
    // Створення дочірнього процесу
    BOOL success = CreateProcessA(
        NULL,                           // Ім'я модуля (NULL - використовувати commandLine)
        const_cast<char*>(commandLine.c_str()), // Командний рядок
        NULL,                           // Атрибути безпеки процесу
        NULL,                           // Атрибути безпеки потоку
        FALSE,                          // Успадкування handle
        0,                              // Прапори створення
        NULL,                           // Блок середовища
        NULL,                           // Поточна директорія
        &si,                            // STARTUPINFO
        &pi                             // PROCESS_INFORMATION
    );
    
    if (!success)
    {
        DWORD error = GetLastError();
        std::cerr << "[Батьківський процес] Помилка створення процесу! Код помилки: " << error << std::endl;
        return 1;
    }
    
    std::cout << "[Батьківський процес] Дочірній процес створено успішно!" << std::endl;
    std::cout << "[Батьківський процес] PID дочірнього процесу: " << pi.dwProcessId << std::endl;
    
    // Закриваємо handle потоку (нам потрібен тільки handle процесу)
    CloseHandle(pi.hThread);
    
    // Крок 2: Опитування стану дочірнього процесу (polling)
    std::cout << "\n[Батьківський процес] Крок 2: Опитування стану дочірнього процесу (polling)..." << std::endl;
    
    DWORD exitCode = STILL_ACTIVE; // STILL_ACTIVE = 259 - процес ще виконується
    int pollCount = 0;
    
    // Цикл опитування: перевіряємо стан процесу, поки він не завершиться
    while (exitCode == STILL_ACTIVE)
    {
        pollCount++;
        
        // Отримуємо код завершення процесу (polling)
        // Якщо процес ще виконується, GetExitCodeProcess поверне STILL_ACTIVE
        if (!GetExitCodeProcess(pi.hProcess, &exitCode))
        {
            DWORD error = GetLastError();
            std::cerr << "[Батьківський процес] Помилка отримання коду завершення! Код помилки: " << error << std::endl;
            CloseHandle(pi.hProcess);
            return 1;
        }
        
        if (exitCode == STILL_ACTIVE)
        {
            std::cout << "[Батьківський процес] Перевірка #" << pollCount << ": процес ще виконується..." << std::endl;
            // Затримка перед наступною перевіркою (polling interval)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    
    std::cout << "[Батьківський процес] Процес завершено після " << pollCount << " перевірок." << std::endl;
    
    // Крок 3: Отримання коду завершення
    std::cout << "\n[Батьківський процес] Крок 3: Отримання коду завершення..." << std::endl;
    std::cout << "[Батьківський процес] Код завершення дочірнього процесу: " << exitCode << std::endl;
    
    // Крок 4: Закриття ресурсів
    std::cout << "\n[Батьківський процес] Крок 4: Закриття ресурсів..." << std::endl;
    CloseHandle(pi.hProcess);
    std::cout << "[Батьківський процес] Handle процесу закрито." << std::endl;
    
    // Крок 5: Фінальний вивід результату
    std::cout << "\n=== РЕЗУЛЬТАТ КОМУНІКАЦІЇ ===" << std::endl;
    std::cout << "Код завершення дочірнього процесу: " << exitCode << std::endl;
    std::cout << "Кількість перевірок (polling): " << pollCount << std::endl;
    std::cout << "=================================" << std::endl;
    
    return 0;
}
