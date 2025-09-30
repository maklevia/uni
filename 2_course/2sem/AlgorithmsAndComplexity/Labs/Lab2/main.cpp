#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <queue>
#include <random>

struct Complex {
    int real, imag;
    Complex(int r = 0, int i = 0) : real(r), imag(i) {}

    double modulus() const {
        return std::sqrt(real * real + imag * imag);
    }

    bool operator<(const Complex& other) const {
        double mod1 = modulus(), mod2 = other.modulus();
        if (mod1 != mod2) return mod1 < mod2;
        return real < other.real;
    }

    bool operator==(const Complex& other) const {
        return real == other.real && imag == other.imag;
    }

    friend std::ostream& operator<<(std::ostream& os, const Complex& c) {
        os << "(" << c.real;
        if (c.imag >= 0) os << " + i*" << c.imag;
        else os << " - i*" << -c.imag;
        os << ")";
        return os;
    }
};

struct Node {
    Complex key;
    int priority;
    Node *left, *right;

    Node(Complex k, int p) : key(k), priority(p), left(nullptr), right(nullptr) {}
};

struct Treap {
    Node* root = nullptr;
    std::set<int> used_priorities;
    std::mt19937 rng{std::random_device{}()};

    int getUniquePriority() {
        std::uniform_int_distribution<int> dist(1, 100);
        int p;
        do { p = dist(rng); } while (used_priorities.count(p));
        used_priorities.insert(p);
        return p;
    }

    void split(Node* t, Complex key, Node*& left, Node*& right) {
        if (!t) {
            left = right = nullptr;
        } else if (t->key < key) {
            split(t->right, key, t->right, right);
            left = t;
        } else {
            split(t->left, key, left, t->left);
            right = t;
        }
    }

    Node* merge(Node* left, Node* right) {
        if (!left || !right) return left ? left : right;
        if (left->priority > right->priority) {
            left->right = merge(left->right, right);
            return left;
        } else {
            right->left = merge(left, right->left);
            return right;
        }
    }

    void insert(Complex key) {
        int priority = getUniquePriority();
        Node* newNode = new Node(key, priority);
        Node *left, *right;
        split(root, key, left, right);
        root = merge(merge(left, newNode), right);
    }

    Node* erase(Node* root, Complex key) {
        if (!root) return nullptr;

        if (root->key == key) {
            Node* newRoot = merge(root->left, root->right);
            delete root;
            return newRoot;
        }

        if (key < root->key)
            root->left = erase(root->left, key);
        else
            root->right = erase(root->right, key);

        return root;
    }

    void erase(Complex key) {
        root = erase(root, key);
    }



    int findLevel(Node* node, Complex key, int level) {
        if (!node) return -1;
        if (node->key == key) return level;
        if (key < node->key) return findLevel(node->left, key, level + 1);
        return findLevel(node->right, key, level + 1);
    }

    void printLevels() {
        if (!root) return;
        std::queue<Node*> q;
        q.push(root);
        int level = 1;
        while (!q.empty()) {
            int size = q.size();
            std::cout << "Level " << level++ << ": ";
            while (size--) {
                Node* node = q.front(); q.pop();
                std::cout << node->key << " [" << node->priority << "] ";
                if (node->left) q.push(node->left);
                if (node->right) q.push(node->right);
            }
            std::cout << "\n";
        }
    }
};

int main() {
    Treap treap;
    int n;
    std::cout << "Enter number of complex numbers: ";
    std::cin >> n;
    std::vector<Complex> values(n);

    std::cout << "Enter " << n << " complex numbers (real imag):\n";
    for (int i = 0; i < n; i++) {
        int real, imag;
        std::cin >> real >> imag;
        values[i] = Complex(real, imag);
    }

    for (const auto& c : values) {
        treap.insert(c);
    }

    while (true) {
        std::cout << "Treap levels:\n";
        treap.printLevels();

        std::cout << "Choose an option: \n";
        std::cout << "1. Insert element\n";
        std::cout << "2. Delete element\n";
        std::cout << "3. Find level of element\n";
        std::cout << "4. Exit\n";
        int choice;
        std::cin >> choice;

        if (choice == 4) break;

        int real, imag;
        switch (choice) {
            case 1:
                std::cout << "Enter complex number to insert (real imag): ";
                std::cin >> real >> imag;
                treap.insert(Complex(real, imag));
                break;
            case 2:
                std::cout << "Enter complex number to delete (real imag): ";
                std::cin >> real >> imag;
                treap.erase(Complex(real, imag));
                break;
            case 3:
                std::cout << "Enter complex number to find (real imag): ";
                std::cin >> real >> imag;
                int level;
                level = treap.findLevel(treap.root, Complex(real, imag), 1);
                if (level == -1) std::cout << "Element not found!\n";
                else std::cout << "Element found at level " << level << "\n";
                break;
            default:
                std::cout << "Invalid choice!\n";
        }
    }

    return 0;
}
