#include <iostream>
#include <vector>
#include <numeric>
using namespace std;

class Fraction {
private:
    int numerator;
    int denominator;

    void simplify() {
        int gcdValue = gcd(abs(numerator), abs(denominator));
        numerator /= gcdValue;
        denominator /= gcdValue;

        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }

public:
    Fraction(int num = 0, int denom = 1) : numerator(num), denominator(denom) {
        simplify();
    }

    Fraction operator+(const Fraction& other) const {
        return Fraction(numerator * other.denominator + other.numerator * denominator,
                        denominator * other.denominator);
    }

    Fraction operator-(const Fraction& other) const {
        return Fraction(numerator * other.denominator - other.numerator * denominator,
                        denominator * other.denominator);
    }

    Fraction operator*(const Fraction& other) const {
        return Fraction(numerator * other.numerator, denominator * other.denominator);
    }

    Fraction operator/(const Fraction& other) const {
        return Fraction(numerator * other.denominator, denominator * other.numerator);
    }

    friend istream& operator>>(istream& is, Fraction& frac) {
        char slash;
        is >> frac.numerator >> slash >> frac.denominator;
        frac.simplify();
        return is;
    }

    friend ostream& operator<<(ostream& os, const Fraction& frac) {
        os << frac.numerator << "/" << frac.denominator;
        return os;
    }
};

class Matrix {
private:
    vector<vector<Fraction>> mat;
    int size;

public:
    Matrix(int n) : size(n), mat(n, vector<Fraction>(n)) {}

    void input() {
        cout << "Enter elements of the matrix in the format `numerator/denominator`:\n";
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cin >> mat[i][j];
            }
        }
    }

    Matrix inverse() const {
        vector<vector<Fraction>> L(size, vector<Fraction>(size));
        vector<vector<Fraction>> U(size, vector<Fraction>(size));
        vector<vector<Fraction>> inv(size, vector<Fraction>(size));

        // LU decomposition
        for (int i = 0; i < size; ++i) {
            for (int j = i; j < size; ++j) {
                Fraction sum(0);
                for (int k = 0; k < i; ++k) {
                    sum = sum + L[i][k] * U[k][j];
                }
                U[i][j] = mat[i][j] - sum;
            }
            for (int j = i; j < size; ++j) {
                if (i == j) {
                    L[i][i] = Fraction(1);
                } else {
                    Fraction sum(0);
                    for (int k = 0; k < i; ++k) {
                        sum = sum + L[j][k] * U[k][i];
                    }
                    L[j][i] = (mat[j][i] - sum) / U[i][i];
                }
            }
        }

        // Find inverse
        for (int i = 0; i < size; ++i) {
            vector<Fraction> y(size);
            for (int j = 0; j < size; ++j) {
                Fraction sum(0);
                for (int k = 0; k < j; ++k) {
                    sum = sum + L[j][k] * y[k];
                }
                y[j] = Fraction((i == j) ? 1 : 0) - sum;
            }
            for (int j = size - 1; j >= 0; --j) {
                Fraction sum(0);
                for (int k = j + 1; k < size; ++k) {
                    sum = sum + U[j][k] * inv[k][i];
                }
                inv[j][i] = (y[j] - sum) / U[j][j];
            }
        }

        return Matrix(inv);
    }

    void display() const {
        for (const auto& row : mat) {
            for (const auto& elem : row) {
                cout << elem << " ";
            }
            cout << endl;
        }
    }

private:
    Matrix(const vector<vector<Fraction>>& data) : mat(data), size(data.size()) {}
};

int main() {
    int n;
    cout << "Enter the size of the matrix (n x n): ";
    cin >> n;

    Matrix A(n);
    A.input();

    cout << "\nThe inverse of matrix A is:\n";
    Matrix invA = A.inverse();
    invA.display();

    return 0;
}


