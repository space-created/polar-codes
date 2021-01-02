#include "GaloisFieldPolynomial.h"
#include "PolarCode.h"

void PolarCode::build_constraint_matrix() {
    vector<vector<u8> > check_matrix = build_bch_check_matrix();

    vector<vector<u8> > f_matrix = kronecker_product();

    vector<vector<u8> > a_matrix = get_a_matrix(f_matrix);

    vector<vector<u8> > a_transposed_matrix = transpose_matrix(a_matrix);

    vector<vector<u8> > h_a_t_product = get_matrix_product(check_matrix, a_transposed_matrix);


    left_triangulation(h_a_t_product);
    triangulation(h_a_t_product);

    for (size_t i = 0; i < h_a_t_product.size() / 2; ++i) {
        swap(h_a_t_product.at(i), h_a_t_product.at(h_a_t_product.size() - i - 1));
    }
    remove_zero_rows(h_a_t_product);

    for (auto & i : constraint_matrix) {
        for (size_t j = 0; j < constraint_matrix.at(0).size(); ++j) {
            cout << (int) i.at(j) << ' ';
        }
        cout << endl;
    }
}

vector<vector<u8> > PolarCode::build_bch_check_matrix() {
    vector<vector<u16> > h(bch_distance - 1, vector<u16>(word_length));
    u16 b = 0;
    for (size_t i = 0; i < bch_distance - 1; ++i) {
        h.at(i).at(0) = 0;
        for (size_t j = 1; j < word_length; ++j) {
            auto gf_pol = GaloisFieldPolynomial(j, n);
            auto mod_poly = GaloisFieldPolynomial(poly, poly.size());
            h.at(i).at(j) = GaloisFieldPolynomial::power_by_modulo(b + i, gf_pol,
                                                             mod_poly).get_representation(); // ^ (b + i)).poly()
        }
    }

    vector<vector<u8> > check_matrix(n * (bch_distance - 1), vector<u8>(word_length, 0));


    for (size_t i = 0; i < bch_distance - 1; ++i) {
        for (size_t j = 0; j < word_length; ++j) {
            u16 temp = h.at(i).at(j);
            for (size_t z = n; z > 0; z--) {
                check_matrix.at(i * n + z - 1).at(j) = temp % 2;
                temp /= 2;
            }
        }
    }
    check_matrix.at(n - 1).at(0) = 1;
    return check_matrix;
}

vector<vector<u8> > PolarCode::get_matrix_product(vector<vector<u8> > &a, vector<vector<u8> > &b) {
    vector<vector<u8> > product(a.size(), vector<u8>(b.at(0).size(), 0));

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.at(0).size(); ++j) {
            for (size_t k = 0; k < a.at(0).size(); ++k) {
                product[i][j] = (product[i][j] + (a[i][k] * b[k][j])) % 2;
            }
        }
    }

    return product;
}

vector<vector<u8> > PolarCode::transpose_matrix(vector<vector<u8> > &matrix) {
    vector<vector<u8> > transposed_matrix(matrix.size(), vector<u8>(matrix[0].size(), 0));
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            transposed_matrix[i][j] = matrix[j][i];
        }
    }
    return transposed_matrix;
}

vector<vector<u8> > PolarCode::get_a_matrix(vector<vector<u8> > &f_matrix) {
    vector<vector<u8> > a_matrix(f_matrix.size(), vector<u8>(f_matrix.size(), 0));
    for (size_t i = 0; i < bit_rev_matrix_order.size(); ++i) {
        a_matrix[i] = move(f_matrix[bit_rev_matrix_order[i]]);
    }
    return a_matrix;
}

vector<vector<u8> > PolarCode::kronecker_product() {
    vector<vector<u8> > arikan_kernel = {{1, 0},
                                         {1, 1}};
    vector<vector<u8> > product = {{1, 0},
                                   {1, 1}};
    for (size_t i = 0; i < n - 1; ++i) {
        vector<vector<u8> > new_product(product.size() * arikan_kernel.size(),
                                        vector<u8>(product[0].size() * arikan_kernel[0].size(), 0));

        for (size_t ia = 0; ia < product.size(); ia++) {
            for (int ja = 0; ja < product[ia].size(); ja++) {
                for (int ib = 0; ib < arikan_kernel.size(); ib++) {
                    for (int jb = 0; jb < arikan_kernel[ib].size(); jb++) {
                        new_product[arikan_kernel.size() * ia + ib][arikan_kernel[ib].size() * ja + jb] =
                                product[ia][ja] * arikan_kernel[ib][jb];
                    }
                }
            }
        }
        product = move(new_product);
    }

    return product;
}

void PolarCode::triangulation(std::vector<std::vector<u8> > &matrix) {
    if (matrix.empty()) {
        return;
    }
    for (int i = matrix.size(); i > 1; i--) {
        unsigned int imax = col_max(matrix, i - 1);
        if ((i - 1) != imax) {
            swap(matrix[i - 1], matrix[imax]);
        }
        for (int j = i - 1; j > 0; j--) {
            u8 mul = matrix[j - 1][i - 1];
            for (int k = i; k > 0; k--) {
                matrix[j - 1][k - 1] = (matrix[j - 1][k - 1] + matrix[i - 1][k - 1] * mul) % 2;
            }
        }
    }
}

void PolarCode::remove_zero_rows(std::vector<std::vector<u8> > &matrix) {
    for (size_t i = 0; i < matrix.size(); ++i) {
        bool is_zero_row = true;
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] == 1) {
                is_zero_row = false;
                break;
            }
        }
        if (!is_zero_row) {
            constraint_matrix.push_back(matrix[i]);
        }
    }

    for (size_t i = 0; i < constraint_matrix.size(); ++i) {
        int pos = -1;
        for (size_t j = 0; j < constraint_matrix[i].size(); ++j) {
            if (pos == -1) {
                if (constraint_matrix[i][j] == 1) {
                    pos = j;
                }
            } else {
                if (constraint_matrix[i][j] == 1) {
                    pos = -1;
                    break;
                }
            }
        }
        if (pos != -1) {
            for (size_t j = 0; j < constraint_matrix.size(); ++j) {
                if (j != i && constraint_matrix[j][pos] == 1) {
                    constraint_matrix[j][pos] = 0;
                }
            }
        }
    }
}

void PolarCode::left_triangulation(std::vector<std::vector<u8> > &matrix) {
    if (matrix.empty())
        return;
    const int num_cols = matrix[0].size();
    for (int i = 0; i < matrix.size() - 1; ++i) {
        unsigned int imax = left_col_max(matrix, i);
        if (i != imax) {
            swap(matrix[i], matrix[imax]);
        }
        for (int j = i + 1; j < matrix.size(); ++j) {
            u8 mul = matrix[j][i];
            for (int k = i; k < num_cols; ++k) {
                matrix[j][k] = (matrix[j][k] + matrix[i][k] * mul) % 2;
            }
        }
    }
}

size_t PolarCode::left_col_max(const std::vector<std::vector<u8> > &matrix, size_t col) {
    u8 max = matrix[col][col];
    int maxPos = col;
    for (int i = col + 1; i < matrix.size(); ++i) {
        u8 element = matrix[i][col];
        if (element > max) {
            max = element;
            maxPos = i;
        }
    }
    return maxPos;
}


size_t PolarCode::col_max(const std::vector<std::vector<u8> > &matrix, size_t col) {
    u8 max = matrix[col][col];
    int maxPos = col;
    for (int i = col; i > 0; i--) {
        u8 element = matrix[i - 1][col];
        if (element > max) {
            max = element;
            maxPos = i - 1;
        }
    }
    return maxPos;
}