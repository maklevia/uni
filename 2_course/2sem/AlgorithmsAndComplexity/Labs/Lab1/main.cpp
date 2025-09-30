#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <cmath>
#include <algorithm>
#include <set>
using namespace std;

int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return abs(a);
}

struct Rational {
    int numerator;
    int denominator;

    Rational(int num, int den) {
        if (den == 0) {
            cerr << "Denominator cannot be zero. Exiting program." << endl;
            exit(0);
        }
        int divisor = gcd(abs(num), abs(den));
        numerator = num / divisor;
        denominator = den / divisor;
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }

    bool operator==(const Rational &other) const {
        return numerator == other.numerator && denominator == other.denominator;
    }

    bool operator<(const Rational &other) const {
        return tie(numerator, denominator) < tie(other.numerator, other.denominator);
    }
};

class SecondaryHashTable {
private:
    vector<vector<Rational>> table;
    int a, b, p, m;

    int hash(const vector<Rational> &values) {
        int sum = 0;
        for (const auto &rat : values) {
            sum += a * rat.numerator + b * rat.denominator;
        }
        return (sum % p) % m;
    }

    void rehash(const vector<vector<Rational>> &values) {
        bool collision;
        do {
            a = 11;
            b = 7;
            p = 101;
            table.assign(m, {});
            collision = false;
            for (const auto &vec : values) {
                int index = hash(vec);
                if (!table[index].empty()) {
                    collision = true;
                    break;
                }
                table[index] = vec;
            }
        } while (collision);
    }

public:
    SecondaryHashTable(const vector<vector<Rational>> &values) {
        m = values.size() * values.size();
        table.resize(m);
        rehash(values);
    }

    bool search(const vector<Rational> &values) {
        int index = hash(values);
        return !table[index].empty() && table[index] == values;
    }

    void display() {
        for (int i = 0; i < m; ++i) {
            cout << "  " << i << ": ";
            if (!table[i].empty()) {
                cout << "(";
                for (const auto &rat : table[i]) {
                    cout << rat.numerator << "/" << rat.denominator << " ";
                }
                cout << ")";
            } else {
                cout << "-";
            }
            cout << endl;
        }
    }
};

class PerfectHashTable {
private:
    vector<vector<vector<Rational>>> primaryTable;
    vector<SecondaryHashTable*> secondaryTables;
    int a, b, p, n;

    int hash(const vector<Rational> &values) {
        int sum = 0;
        for (const auto &rat : values) {
            sum += a * rat.numerator + b * rat.denominator;
        }
        return (sum % p) % n;
    }

    void rehash() {
        int attempts = 0;
        do {
            a = 7;
            b = 11;
            p = 101;
            vector<vector<vector<Rational>>> newTable(n);
            for (const auto &bucket : primaryTable) {
                for (const auto &vec : bucket) {
                    int index = hash(vec);
                    newTable[index].push_back(vec);
                }
            }
            bool valid = true;
            for (const auto &bucket : newTable) {
                if (bucket.size() > 1) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                primaryTable = move(newTable);
                break;
            }
            attempts++;
        } while (attempts < 100);
        secondaryTables.clear();
        for (const auto &bucket : primaryTable) {
            if (bucket.size() > 1) {
                secondaryTables.push_back(new SecondaryHashTable(bucket));
            } else {
                secondaryTables.push_back(nullptr);
            }
        }
    }

public:
    PerfectHashTable(int size) : n(size) {
        primaryTable.resize(n);
        srand(time(0));
        rehash();
    }

    void insert(const vector<Rational> &vec) {
        int index = hash(vec);
        for (const auto &existing : primaryTable[index]) {
            if (existing == vec) return;
        }
        primaryTable[index].push_back(vec);
        rehash();
    }

    bool search(const vector<Rational> &vec) {
        int index = hash(vec);
        for (const auto &existing : primaryTable[index]) {
            if (existing == vec) return true;
        }
        if (secondaryTables[index]) {
            return secondaryTables[index]->search(vec);
        }
        return false;
    }

    void display() {
        for (int i = 0; i < n; ++i) {
            cout << "Index " << i << ": ";
            if (!primaryTable[i].empty() and !secondaryTables[i]) {
                cout << "(";
                for (const auto &rat : primaryTable[i][0]) {
                    cout << rat.numerator << "/" << rat.denominator << " ";
                }
                cout << ") ";
            } else if (!secondaryTables[i]) {
                cout << "_";
            }
            if (secondaryTables[i]) {
                cout << "\nSecondary Hash Table:" << endl;
                secondaryTables[i]->display();
            }
            cout << endl;
        }
    }
};
int main() {
    int count;
    cout << "Enter the number of sets of rational numbers: ";
    cin >> count;

    set<vector<Rational>> uniqueSets;
    for (int i = 0; i < count; ++i) {
        int numElements;
        cout << "Enter number of elements in set " << i + 1 << ": ";
        cin >> numElements;
        vector<Rational> temp;
        for (int j = 0; j < numElements; ++j) {
            int num, den;
            cout << "Enter numerator and denominator for element " << j + 1 << ": ";
            cin >> num >> den;
            temp.emplace_back(num, den);
        }
        uniqueSets.insert(temp);
    }
    PerfectHashTable hashTable(uniqueSets.size());
    for (const auto &set : uniqueSets) {
        hashTable.insert(set);
    }
    hashTable.display();

    vector<Rational> query;
    int numElements;
    cout << "Enter number of elements in the set to search: ";
    cin >> numElements;
    for (int i = 0; i < numElements; ++i) {
        int num, den;
        cout << "Enter numerator and denominator for element " << i + 1 << ": ";
        cin >> num >> den;
        query.emplace_back(num, den);
    }

    if (hashTable.search(query)) {
        cout << "The set exists in the hash table." << endl;
    } else {
        cout << "The set was not found in the hash table." << endl;
    }

    cout << "Enter number of elements in the set to search: ";
    cin >> numElements;
    for (int i = 0; i < numElements; ++i) {
        int num, den;
        cout << "Enter numerator and denominator for element " << i + 1 << ": ";
        cin >> num >> den;
        query.emplace_back(num, den);
    }

    if (hashTable.search(query)) {
        cout << "The set exists in the hash table." << endl;
    } else {
        cout << "The set was not found in the hash table." << endl;
    }
    return 0;
}