#include <iostream>
#include <vector>
#include <string>
#include <cstdint>

const int64_t mod = 1e9 + 7;
const int base = 101;
const int base_row = 103;

// Function to read a grid from the input
void readGrid(std::vector<std::string>& grid, int size, const std::string& name) {
    std::cout << "Enter the " << name << " grid:\n";
    for (int i = 0; i < size; ++i) {
        std::cin >> grid[i];
    }
}

// Function to precompute powers for rolling hash
void precomputePowers(std::vector<int64_t>& pow_base, std::vector<int64_t>& pow_base_row, int max_len) {
    pow_base[0] = pow_base_row[0] = 1;
    for (int i = 1; i < max_len; ++i) {
        pow_base[i] = (pow_base[i - 1] * base) % mod;
        pow_base_row[i] = (pow_base_row[i - 1] * base_row) % mod;
    }
}

// Function to compute the hash of the pattern grid
int64_t computePatternHash(const std::vector<std::string>& pattern, int m) {
    int64_t pattern_hash = 0;
    for (int i = 0; i < m; ++i) {
        int64_t row_hash = 0;
        for (int j = 0; j < m; ++j) {
            row_hash = (row_hash * base + pattern[i][j]) % mod;
        }
        pattern_hash = (pattern_hash * base_row + row_hash) % mod;
    }
    return pattern_hash;
}

// Function to compute row hashes for the text grid
void computeTextRowHashes(const std::vector<std::string>& text, std::vector<std::vector<int64_t>>& text_row_hashes, int n, int m, const std::vector<int64_t>& pow_base) {
    for (int i = 0; i < n; ++i) {
        int64_t h = 0;
        for (int j = 0; j < m; ++j) {
            h = (h * base + text[i][j]) % mod;
        }
        text_row_hashes[i][0] = h;
        for (int j = 1; j <= n - m; ++j) {
            h = ((h - text[i][j - 1] * pow_base[m - 1] % mod + mod) * base + text[i][j + m - 1]) % mod;
            text_row_hashes[i][j] = h;
        }
    }
}

// Function to search for the pattern in the text grid
void searchPattern(const std::vector<std::string>& text, const std::vector<std::string>& pattern, int n, int m, const std::vector<std::vector<int64_t>>& text_row_hashes, int64_t pattern_hash, const std::vector<int64_t>& pow_base_row) {
    bool found = false;

    for (int j = 0; j <= n - m; ++j) {
        int64_t block_hash = 0;
        for (int k = 0; k < m; ++k) {
            block_hash = (block_hash * base_row + text_row_hashes[k][j]) % mod;
        }

        if (block_hash == pattern_hash) {
            bool match = true;
            for (int x = 0; x < m && match; ++x) {
                for (int y = 0; y < m; ++y) {
                    if (text[x][j + y] != pattern[x][y]) {
                        match = false;
                        break;
                    }
                }
            }
            if (match) {
                std::cout << "Pattern found at position (0, " << j << ")\n";
                found = true;
            }
        }

        for (int i = 1; i <= n - m; ++i) {
            block_hash = ((block_hash - text_row_hashes[i - 1][j] * pow_base_row[m - 1] % mod + mod) * base_row
                + text_row_hashes[i + m - 1][j]) % mod;

            if (block_hash == pattern_hash) {
                bool match = true;
                for (int x = 0; x < m && match; ++x) {
                    for (int y = 0; y < m; ++y) {
                        if (text[i + x][j + y] != pattern[x][y]) {
                            match = false;
                            break;
                        }
                    }
                }
                if (match) {
                    std::cout << "Pattern found at position (" << i << ", " << j << ")\n";
                    found = true;
                }
            }
        }
    }

    if (!found) {
        std::cout << "Pattern not found in the text grid.\n";
    }
}

int main() {
    int n, m;

    std::cout << "Enter the size of the text grid (n): ";
    std::cin >> n;
    std::vector<std::string> text(n);
    readGrid(text, n, "text");

    std::cout << "Enter the size of the pattern grid (m): ";
    std::cin >> m;
    std::vector<std::string> pattern(m);
    readGrid(pattern, m, "pattern");

    int max_len = std::max(n, m) + 1;
    std::vector<int64_t> pow_base(max_len), pow_base_row(max_len);
    precomputePowers(pow_base, pow_base_row, max_len);

    int64_t pattern_hash = computePatternHash(pattern, m);

    std::vector<std::vector<int64_t>> text_row_hashes(n, std::vector<int64_t>(n - m + 1));
    computeTextRowHashes(text, text_row_hashes, n, m, pow_base);

    searchPattern(text, pattern, n, m, text_row_hashes, pattern_hash, pow_base_row);

    return 0;
}

