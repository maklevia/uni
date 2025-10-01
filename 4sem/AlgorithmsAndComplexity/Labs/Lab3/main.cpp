#include <iostream>
#include <cmath>
#include <limits>

struct Complex {
    double real, imag;

    Complex(double r = 0, double i = 0) : real(r), imag(i) {}

    double modulus() const {
        return std::sqrt(real * real + imag * imag);
    }

    bool operator<(const Complex& other) const {
        if (modulus() != other.modulus())
            return modulus() < other.modulus();
        return real < other.real;
    }

    bool operator==(const Complex& other) const {
        return real == other.real && imag == other.imag;
    }

    friend std::ostream& operator<<(std::ostream& os, const Complex& c) {
        os << c.real << (c.imag >= 0 ? "+" : "") << c.imag << "i";
        return os;
    }
};

struct BinomialNode {
    Complex key;
    int degree;
    BinomialNode* parent;
    BinomialNode* child;
    BinomialNode* sibling;

    BinomialNode(const Complex& k) : key(k), degree(0), parent(nullptr), child(nullptr), sibling(nullptr) {}
};

class BinomialHeap {
    BinomialNode* head;

public:
    BinomialHeap() : head(nullptr) {}

    void insert(const Complex& key) {
        BinomialHeap temp;
        temp.head = new BinomialNode(key);
        head = merge(head, temp.head);
    }

    Complex findMin() const {
        if (!head) throw std::runtime_error("Heap is empty");
        BinomialNode* y = nullptr;
        BinomialNode* x = head;
        Complex min = {std::numeric_limits<double>::max(), 0};
        while (x) {
            if (x->key < min) {
                min = x->key;
                y = x;
            }
            x = x->sibling;
        }
        return y->key;
    }

    void deleteMin() {
        if (!head) return;

        BinomialNode *prevMin = nullptr, *minNode = head, *prev = nullptr;
        Complex min = head->key;
        BinomialNode* curr = head->sibling;
        prev = head;

        while (curr) {
            if (curr->key < min) {
                min = curr->key;
                prevMin = prev;
                minNode = curr;
            }
            prev = curr;
            curr = curr->sibling;
        }

        if (prevMin) prevMin->sibling = minNode->sibling;
        else head = minNode->sibling;

        BinomialNode* child = minNode->child;
        BinomialNode* rev = nullptr;
        while (child) {
            BinomialNode* next = child->sibling;
            child->sibling = rev;
            child->parent = nullptr;
            rev = child;
            child = next;
        }

        head = merge(head, rev);
        delete minNode;
    }

    void decreaseKey(BinomialNode* x, const Complex& k) {
        if (x->key < k) {
            throw std::runtime_error("New key is bigger than initial key");
        }

        x->key = k;
        BinomialNode* y = x;
        BinomialNode* z = y->parent;

        while (z && y->key < z->key) {
            std::swap(y->key, z->key);
            y = z;
            z = y->parent;
        }
    }


    void print() const {
        printTree(head);
    }

    BinomialNode* getHead() const { return head; }

private:
    static BinomialNode* merge(BinomialNode* h1, BinomialNode* h2) {
        if (!h1) return h2;
        if (!h2) return h1;

        BinomialNode* newHead = nullptr;
        BinomialNode** pos = &newHead;

        while (h1 && h2) {
            if (h1->degree <= h2->degree) {
                *pos = h1;
                h1 = h1->sibling;
            } else {
                *pos = h2;
                h2 = h2->sibling;
            }
            pos = &(*pos)->sibling;
        }
        *pos = h1 ? h1 : h2;

        return unionTrees(newHead);
    }

    static BinomialNode* unionTrees(BinomialNode* head) {
        if (!head) return nullptr;

        BinomialNode *prev = nullptr, *curr = head, *next = head->sibling;

        while (next) {
            if (curr->degree != next->degree || (next->sibling && next->sibling->degree == curr->degree)) {
                prev = curr;
                curr = next;
            } else if (curr->key < next->key) {
                curr->sibling = next->sibling;
                linkTrees(curr, next);
            } else {
                if (!prev) head = next;
                else prev->sibling = next;
                linkTrees(next, curr);
                curr = next;
            }
            next = curr->sibling;
        }

        return head;
    }

    static void linkTrees(BinomialNode* root, BinomialNode* child) {
        child->parent = root;
        child->sibling = root->child;
        root->child = child;
        root->degree++;
    }

    void printTree(BinomialNode* node) const {
        while (node) {
            std::cout << "B" << node->degree << ":\n";
            printSubtree(node, 1);
            node = node->sibling;
        }
    }

    void printSubtree(BinomialNode* node, int depth) const {
        for (int i = 0; i < depth - 1; ++i) std::cout << "  ";
        std::cout << node->key << "\n";
        BinomialNode* child = node->child;
        while (child) {
            printSubtree(child, depth + 1);
            child = child->sibling;
        }
    }
};

int main() {
    BinomialHeap heap;
    int n;
    std::cout << "Enter number of complex numbers: ";
    std::cin >> n;
    std::cout << "Enter complex numbers (real imag):\n";
    for (int i = 0; i < n; ++i) {
        double r, im;
        std::cin >> r >> im;
        heap.insert(Complex(r, im));
    }
    std::cout << "\nInitial heap:\n";
    heap.print();

    int choice;
    do {
        std::cout << "\nMenu:\n1. Insert\n2. Delete Min\n3. Decrease Key\n0. Exit\nChoice: ";
        std::cin >> choice;

        if (choice == 1) {
            double r, im;
            std::cout << "Enter complex number to insert: ";
            std::cin >> r >> im;
            heap.insert(Complex(r, im));
            heap.print();
        } else if (choice == 2) {
            heap.deleteMin();
            heap.print();
        } else if (choice == 3) {
            std::cout << "Enter current complex number (real imag): ";
            double currR, currI;
            std::cin >> currR >> currI;
            Complex currentKey(currR, currI);

            std::cout << "Enter new decreased value (real imag): ";
            double newR, newI;
            std::cin >> newR >> newI;
            Complex newKey(newR, newI);

            BinomialNode* nodeToDecrease = nullptr;

            std::function<void(BinomialNode*)> findNode = [&](BinomialNode* curr) {
                if (!curr || nodeToDecrease) return;
                if (curr->key == currentKey) {
                    nodeToDecrease = curr;
                    return;
                }
                findNode(curr->child);
                findNode(curr->sibling);
            };

            findNode(heap.getHead());

            if (nodeToDecrease) {
                try {
                    heap.decreaseKey(nodeToDecrease, newKey);
                    std::cout << "Key decreased successfully.\n";
                    heap.print();
                } catch (const std::runtime_error& e) {
                    std::cout << e.what() << "\n";
                }
            } else {
                std::cout << "Node with the given current value not found.\n";
            }
        } else if (choice == 0) {
            break;
        }

    } while (choice != 0);

    return 0;
}

