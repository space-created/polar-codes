#include "GaloisFieldPolynomial.h"
#include "PolarCode.h"
/*
 * Method for building EBCH code constraint matrix (https://arxiv.org/abs/1511.01646 p.4)
 */
void PolarCode::build_constraint_matrix() {
    vector<vector<u8> > check_matrix = build_bch_check_matrix();

    vector<vector<u8> > f_matrix = kronecker_product();

    vector<vector<u8> > a_matrix = get_a_matrix(f_matrix);

    vector<vector<u8> > a_transposed_matrix = transpose_matrix(a_matrix);
    vector<vector<u8> > h_a_t_product = get_matrix_product(check_matrix, a_transposed_matrix);

    right_triangulation(h_a_t_product);


    sort_by_right_one(h_a_t_product);
    matrix_reduction(h_a_t_product);
    remove_zero_rows(h_a_t_product);


    sort_by_right_one(constraint_matrix);
    cout << "First\n";
    cout << (int) m << ' ' << constraint_matrix.size() << ' ' << constraint_matrix.at(0).size()  << '\n';
    cout << '\n';
    for (auto & i : constraint_matrix) {
        for (size_t j = 0; j < constraint_matrix.at(0).size(); ++j) {
            cout << (int) i.at(j) << ' ';
        }
        cout << endl;
    }
    J.resize(constraint_matrix.size(), -1);
    for (size_t i = 0; i < J.size(); ++i) {
        J.at(i) = find_the_most_right_one_pos(constraint_matrix.at(i));
    }
    initialize_frozen_bits();
    insert_least_reliable_rows_for_freezing_bits();
    sort_by_right_one(constraint_matrix);
    cout << (int) m << ' ' << constraint_matrix.size() << ' ' << constraint_matrix.at(0).size() << ' ' << J.at(0) << '\n';
    cout << '\n';
    for (auto & i : constraint_matrix) {
        for (size_t j = 0; j < constraint_matrix.at(0).size(); ++j) {
            cout << (int) i.at(j) << ' ';
        }
        cout << endl;
    }

    if (shoud_read_matrix) {
        for (int i = 0; i < constraint_matrix.size(); ++i) {
            int pos = find_the_most_right_one_pos(constraint_matrix.at(i));
            for (int j = 0; j < pos; ++j) {
                    constraint_matrix.at(i).at(j) = 0;
            }
        }
        int s;
        cin >> s;
//        cout << s << flush;
        for (int i = 0; i < s; ++i) {
            int size_c, frozen_b;
            cin >> size_c >> frozen_b;
            int m_i = 0;
            for (int ii = 0; ii < constraint_matrix.size(); ++ii) {
                if (constraint_matrix.at(ii).at(frozen_b) == 1) {
                    m_i = ii;
                    break;
                }
            }
            for (int ii = 0; ii < constraint_matrix.at(m_i).size(); ++ii) {
                if (ii != frozen_b) {
                    constraint_matrix.at(m_i).at(ii) = 0;
                }
            }
            for (int j = 0; j < size_c; ++j) {
                int ind = 0;
                cin >> ind;
                constraint_matrix.at(m_i).at(ind) = 1;
            }
        }
    }


    J.resize(constraint_matrix.size(), -1);
    // T - is the map which maps the position of the most right 1 to row number of V: e.g. 960 -> 0
    T_arr.resize(constraint_matrix[0].size(), -1);
//    cout << "\nJ: ";
    for (size_t i = 0; i < J.size(); ++i) {

        int one_pos = find_the_most_right_one_pos(constraint_matrix.at(i));
        J.at(i) = one_pos;
        T_arr.at(one_pos) = i;
//        cout << J.at(i) << ' ';
    }
//    cout << "\n";
}

void PolarCode::insert_least_reliable_rows_for_freezing_bits() {
    int amount_to_freeze = static_frozen_channels.size();
    for (int i = 0; i < static_frozen_channels.size(); ++i) {
        vector<u8> V_row_to_insert(constraint_matrix.at(0).size(), 0);
        V_row_to_insert.at(static_frozen_channels.at(i)) = 1;
        if (amount_to_freeze <= q) {
            for (int j = 0; j < static_frozen_channels.at(i); ++j) {
                if (frozen_bits.at(j) == 0) {
                    V_row_to_insert.at(j) = rand() % 2;
                }
            }
        }
        constraint_matrix.push_back(V_row_to_insert);
        amount_to_freeze--;
    }
}

vector<vector<u8> > PolarCode::build_bch_check_matrix() {
    vector<vector<u16> > h(bch_distance - 1, vector<u16>(word_length));
    u16 b = 0;
    for (size_t i = 0; i < bch_distance - 1; ++i) {
        h.at(i).at(0) = 0;
        for (size_t j = 1; j < word_length; ++j) {
            auto gf_pol = GaloisFieldPolynomial(j, m);
            auto mod_poly = GaloisFieldPolynomial(poly, poly.size());
            h.at(i).at(j) = GaloisFieldPolynomial::power_by_modulo(b + i, gf_pol,
                                                             mod_poly).get_representation(); // ^ (b + i)).poly()
        }
    }

    vector<vector<u8> > check_matrix(m * (bch_distance - 1), vector<u8>(word_length, 0));

    for (size_t i = 0; i < bch_distance - 1; ++i) {
        for (size_t j = 0; j < word_length; ++j) {
            u16 temp = h.at(i).at(j);
            for (size_t z = m; z > 0; z--) {
                check_matrix.at(i * m + z - 1).at(j) = temp % 2;
                temp /= 2;
            }
        }
    }
    check_matrix.at(m - 1).at(0) = 1;
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
    for (size_t i = 0; i < m - 1; ++i) {
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
}


size_t PolarCode::right_col_max(const std::vector<std::vector<u8> > &matrix, size_t row, size_t col) {
    int max = matrix[row][col];
    int maxPos = row;
    for (int i = row + 1; i < matrix.size(); ++i) {
        int element = matrix[i][col];
        if (element > max) {
            max = element;
            maxPos = i;
        }
    }
    return maxPos;
}

int PolarCode::count_right_zero_columns(std::vector<std::vector<u8> > &matrix, int row_pos, int pos) {
    int count = 0;
    for (int i = pos; i >= 0; --i) {
        for (int j = row_pos; j < matrix.size(); ++j) {
            if (matrix[j][i] == 1) {
                return count;
            }
        }
        count++;
    }
    return count;
}

void PolarCode::right_triangulation(std::vector<std::vector<u8> > &matrix) {
    if (matrix.empty())
        return;
    int num_col = matrix[0].size();
    for (int i = 0; i < matrix.size() - 1; ++i) {
        num_col -= count_right_zero_columns(matrix, i, num_col - 1 - i);
        int imax = right_col_max(matrix, i, num_col - 1 - i);
        if (i != imax) {
            swap(matrix[i], matrix[imax]);
        }
        for (int j = i + 1; j < matrix.size(); ++j) {
            int mul = matrix[j][num_col - 1 - i];
            for (int k = num_col - 1 - i; k >= 0; k--) {
                matrix[j][k] = (matrix[j][k] + matrix[i][k] * mul) % 2;
            }
        }
    }
}

int PolarCode::find_the_most_right_one_pos(std::vector<u8> &row_vector) {
    int pos = row_vector.size() - 1;
    while (pos > -1) {
        if (row_vector[pos]) {
            return pos;
        }
        --pos;
    }
    return pos;
}

void PolarCode::matrix_reduction(vector<vector<u8>>& matrix) {
    for (int i = matrix.size() - 1; i >= 1; --i) {
        int right_one_pos = find_the_most_right_one_pos(matrix[i]);
        if (right_one_pos != -1) {
            for (int j = i - 1; j >= 0; --j) {
                int count_zero = 0;
                int count_one = 0;
                for (int k = 0; k <= right_one_pos; ++k) {
                    if (matrix[j][k] == matrix[i][k]) {
                        if (matrix[i][k] == 1) {
                            count_zero++;
                        }
                    } else {
                        if (matrix[i][k] == 1) {
                            count_one++;
                        }
                    }
                }
                if (count_zero > count_one) {
                    for (int k = 0; k < matrix[0].size(); ++k) {
                        matrix[j][k] = (matrix[j][k] + matrix[i][k]) % 2;
                    }
                }
            }
        }
    }
}

void PolarCode::sort_by_right_one(vector<vector<u8>>& matrix) {
    std::sort(matrix.begin(), matrix.end(), [this](vector<u8> & a, vector<u8> & b)
    {
        return (find_the_most_right_one_pos(a) > find_the_most_right_one_pos(b));
    });
}