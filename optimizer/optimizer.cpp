// Include the header file required to use the generic
// interface for the C++ shared library generated by the
// MATLAB Compiler SDK.
#include "MatlabCppSharedLib.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <chrono>
#include <algorithm>
#include <random>
#include <set>

#define u8 uint8_t
#define u16 uint16_t
#define ll long long

using namespace std;
using namespace std::chrono;

namespace mc = matlab::cpplib;
namespace md = matlab::data;
std::shared_ptr<mc::MATLABApplication> setup()
{
	auto mode = mc::MATLABApplicationMode::IN_PROCESS;
	// Specify MATLAB startup options
	std::vector<std::u16string> options = {};
	std::shared_ptr<mc::MATLABApplication> matlabApplication = mc::initMATLABApplication(mode, options);
	return matlabApplication;
}

ll calcACount = 0;
ll computeWEFCountNaive = 0;
ll computeWEFCount = 0;

typedef complex<double> cd;
const double PI = acos(-1);

int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

void fft(vector<cd> & a, bool invert) {
    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;

    for (int i = 0; i < n; i++) {
        if (i < reverse(i, lg_n))
            swap(a[i], a[reverse(i, lg_n)]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd & x : a)
            x /= n;
    }
}

vector<ll> multiply(vector<ll> const& a, vector<ll> const& b) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    int final_size = a.size() + b.size();
    while (n < final_size) {
        n <<= 1;
    }
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }
    fft(fa, true);

    vector<ll> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = round(fa[i].real());
    }
    vector<ll> res = vector<ll>(result.begin(), result.end() - (result.size() - final_size + 1));
    return res;
}

vector<ll> multiply2(vector<ll> const& a, vector<ll> const& b) {
    vector<ll> prod(a.size() + b.size() - 1);
    for (int i = 0; i< prod.size(); i++)
        prod[i] = 0;

    for (int i=0; i<a.size(); i++) {
        for (int j=0; j<b.size(); j++)
            prod[i+j] += a[i] * b[j];
    }

    return prod;
}

vector<ll> add(const vector<ll>& a, const vector<ll> & b) {
    vector<ll> sum(max(a.size(), b.size()), 0);
    for (int i = 0; i < sum.size(); ++i) {
        if (i < a.size() && i < b.size()) {
            sum.at(i) = a.at(i) + b.at(i);
        } else if (i >= a.size()) {
            sum.at(i) = b.at(i);
        } else if (i >= b.size()) {
            sum.at(i) = a.at(i);
        }
    }
    return sum;
}

//vector<ll> multiply_constant(ll constant, vector<ll> const& a) {
//    vector<ll> result(a.size(), 0);
//    for (int i = 0; i < a.size(); ++i) {
//        result.at(i) = constant * a.at(i);
//    }
//    return result;
//}


pair <vector<u8>, vector<u8>> get_even_odd_vectors(vector <u8> u) {
    vector <u8> u_even(u.size() / 2, 0);
    vector <u8> u_odd(u.size() / 2, 0);
    for (int i = 0; i < u.size(); ++i) {
        if (i % 2 == 0) {
            u_even.at(i / 2) = u.at(i);
        } else {
            u_odd.at(i / 2) = u.at(i);
        }
    }
    return {u_even, u_odd};
}

pair<vector<ll>, vector<ll> > calcA(int n, vector <u8> u) {
    calcACount++;
    if (n == 1) {
        return {{1, 0},
                {0, 1}};
    }
    int i = u.size();
    if (i % 2 == 0) {
        pair <vector<u8>, vector<u8>> even_odd_u = get_even_odd_vectors(u);
        for (int i = 0; i < even_odd_u.first.size(); ++i) {
            even_odd_u.first.at(i) = even_odd_u.first.at(i) ^ even_odd_u.second.at(i);
        }
        pair<vector<ll>, vector<ll> > f = calcA(n / 2, even_odd_u.first);
        pair<vector<ll>, vector<ll> > g = calcA(n / 2, even_odd_u.second);


        return {add(multiply2(f.first, g.first),
                    multiply2(f.second, g.second)),
                add(multiply2(f.first, g.second),
                    multiply2(f.second, g.first))};
    } else {
        vector<u8> u_short = vector<u8>(u.begin(), u.begin() + i - 1);
        pair <vector<u8>, vector<u8>> even_odd_u = get_even_odd_vectors(u_short);
        for (int i = 0; i < even_odd_u.first.size(); ++i) {
            even_odd_u.first.at(i) = even_odd_u.first.at(i) ^ even_odd_u.second.at(i);
        }
        pair<vector<ll>, vector<ll> > f = calcA(n / 2, even_odd_u.first);
        pair<vector<ll>, vector<ll> > g =  calcA(n / 2, even_odd_u.second);
        if (u.at(u.size() - 1) == 0) {
            return {multiply2(f.first, g.first),
                    multiply2(f.second, g.second)};
        } else {
            return {multiply2(f.first, g.second),
                    multiply2(f.second, g.first)};
        }
    }
}

//void traverseWithRedBitsAndDynamicConstraints(vector<pair<bool, vector<int> > >& dynamic_constraints,
//                                              int n, int *arr, int i,
//        vector<ll>& a, int last_frozen_pos, int word_size, vector<int> red_indexes) {
//    if (i == n) {
//        vector<int> temp(n, 0);
//        for (int j = 0; j < n; ++j) {
//            temp.at(j) = arr[j];
//        }
//        vector<u8> u_short = vector<u8>(temp.begin(), temp.begin() + last_frozen_pos);
//        pair <vector<ll>, vector<ll>> f = calcA(word_size, u_short);
//        if (temp.at(last_frozen_pos) == 0) {
//            a = add(a, f.first);
//        } else {
//            a = add(a, f.second);
//        }
//        return;
//    }
//    if (!dynamic_constraints.at(i).first) {
//        if (red_indexes.at(i) == -1) {
//            arr[i] = 0;
//            traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, n, arr,
//                    i + 1, a, last_frozen_pos, word_size, red_indexes);
//
//            arr[i] = 1;
//            traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, n, arr,
//                    i + 1, a, last_frozen_pos, word_size, red_indexes);
//        } else {
//            arr[i] = red_indexes.at(i);
//            traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, n, arr,
//                    i + 1, a, last_frozen_pos, word_size, red_indexes);
//        }
//    } else {
//        arr[i] = 0;
//        if (dynamic_constraints.at(i).second.empty()) {
//            traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, n, arr, i + 1,
//                                                     a, last_frozen_pos, word_size, red_indexes);
//        } else {
//            for (int ii = 0; ii < dynamic_constraints.at(i).second.size(); ++ii) {
//                arr[i] = (arr[i] + arr[dynamic_constraints.at(i).second.at(ii)]) % 2;
//            }
//            traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, n, arr, i + 1,
//                                                     a, last_frozen_pos, word_size, red_indexes);
//        }
//    }
//}

vector<ll> computeWEFNaive(int n, int last_frozen_pos,
                           vector<pair<bool, vector<int> > >& dynamic_constraints, const vector<int>& red_indexes_values) {
    computeWEFCountNaive++;

    vector<ll> a(2 * n + 1, 0);
    //    int arr[last_frozen_pos + 1];
    //    traverseWithRedBitsAndDynamicConstraints(dynamic_constraints, last_frozen_pos + 1, arr,
    //                                             0, a, last_frozen_pos, n, red_indexes_values);
    ll cur_number_red_bits = 0;
    for (int i = 0; i < red_indexes_values.size(); ++i) {
        if (red_indexes_values.at(i) == -1) {
            cur_number_red_bits++;
        }
    }
    for (ll i = 0; i < (pow(2, cur_number_red_bits)); ++i) {
        vector<int> cur_vec(n, 0);
        ll temp = i;

        for (int j = 0; j < n; ++j) {
            if (red_indexes_values.at(j) == -1) {
                cur_vec.at(j) = temp % 2;
                temp /= 2;
            } else if (red_indexes_values.at(j) != -2) {
                cur_vec.at(j) = red_indexes_values.at(j);
            }
        }
        for (int j = 0; j < n; ++j) {
            if (!dynamic_constraints.at(j).second.empty()) {
                for (int ii = 0; ii < dynamic_constraints.at(j).second.size(); ++ii) {
                    cur_vec.at(j) = (cur_vec.at(j) + cur_vec.at(dynamic_constraints.at(j).second.at(ii))) % 2;
                }
            }
        }
        vector<u8> u_short = vector<u8>(cur_vec.begin(), cur_vec.begin() + last_frozen_pos);
        pair<vector<ll>, vector<ll> > f = calcA(n, u_short);

        if (cur_vec.at(last_frozen_pos) == 0) {
            a = add(a, f.first);
        } else {
            a = add(a, f.second);
        }
    }

    return a;
}

int compareTwoMonomials(vector<int> monomial_a, vector<int> monomial_b) {
    int var_num_a = 0;
    int var_num_b = 0;
    for (int i = 0; i < monomial_a.size(); ++i) {
        if (monomial_a.at(i) == 1) {
            var_num_a++;
        }
        if (monomial_b.at(i) == 1) {
            var_num_b++;
        }
    }
    if (abs(var_num_a - var_num_b) > 1) {
        return 0;
    } else if (abs(var_num_a - var_num_b) == 1) {
        int flag = 0;
        int res = 1;
        for (int i = monomial_a.size() - 1; i >= 0; --i) {
            if (monomial_a.at(i) > monomial_b.at(i)) {
                if (flag == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = 1;
                }
                flag++;
            } else if (monomial_a.at(i) < monomial_b.at(i)) {
                if (flag == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = -1;
                }
                flag++;
            }
        }
        return res;
    } else {
        int flag = 0;
        int res = 1;
        for (int i = monomial_a.size() - 1; i >= 0; --i) {
            if (monomial_a.at(i) > monomial_b.at(i)) {
                if (flag == 2) {
                    return 0;
                }
                if (flag == 1 && res == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = 1;
                }
                flag++;
            } else if (monomial_a.at(i) < monomial_b.at(i)) {
                if (flag == 2) {
                    return 0;
                }
                if (flag == 1 && res == -1) {
                    return 0;
                }
                if (flag == 0) {
                    res = -1;
                }
                flag++;
            }
        }
        return res;
    }
}

vector<vector<int>> get_monomials_order(vector<vector<int>> monomials, vector<int> info_indexes) {
    vector<vector<int>> monomials_order(monomials.size(), vector<int>(monomials.size(), 0));
    for (int i = 0; i < monomials.size(); ++i) {
        for (int j = 0; j < monomials.size(); ++j) {
            if (i == j) {
                monomials_order.at(i).at(j) = 0;
            } else {
                monomials_order.at(i).at(j) = compareTwoMonomials(monomials.at(i), monomials.at(j));
            }
        }
    }
    ll lambda = 0;
    for (int ii = 0; ii < info_indexes.size(); ++ii) {
        int ind = info_indexes.at(ii);
        int ff = 0;
        for (int i = ind + 1; i < monomials_order.at(ind).size(); ++i) {
            if (monomials_order.at(ind).at(i) == -1) {

                ff++;
            }
        }
        lambda += pow(2, ff);
//        cout << ind << " " << ff << "\n";
    }
//    cout << "LAMBDA= " << lambda << "\n";
    return monomials_order;
}

vector<ll> computeWEF(int n, int s,
                      vector<vector<int>>& monomials_order,
                      vector<pair<int, bool>>& red_indexes,
                      vector<pair<bool, vector<int> > >& dynamic_constraints, vector<int> red_indexes_values) {
    computeWEFCount++;
    bool flag = false;
    for (int i = 0; i < red_indexes.size(); ++i) {
        flag = flag | red_indexes.at(i).second;
    }
    if (!flag) {
        vector<u8> u_zero(s, 0);
        pair<vector<ll>, vector<ll> > f = calcA(n, u_zero);
        return f.first;
    } else {
        int f, ii;
        for (int i = 0; i < red_indexes.size(); ++i) {
            if (red_indexes.at(i).second) {
                f = red_indexes.at(i).first;
                ii = i;
                break;
            }
        }
        vector<int> s_set;

        for (int i = 0; i < red_indexes.size(); ++i) {
            if (red_indexes.at(i).second && monomials_order.at(f).at(red_indexes.at(i).first) == 1) {
                s_set.push_back(red_indexes.at(i).first);
            }
        }

//        vector<int> new_red_indexes;
//        for (int i = 0; i < red_indexes.size(); ++i) {
//            if (red_indexes.at(i) != f) {
//                new_red_indexes.push_back(red_indexes.at(i));
//            }
//        }
        red_indexes.at(ii).second = 0;

        red_indexes_values.at(f) = 0;
        vector<ll> result_b = computeWEF(n, s,
                                         monomials_order, red_indexes, dynamic_constraints,
                                         red_indexes_values);
        red_indexes_values.at(f) = 1;

        for (int i = 0; i < s_set.size(); ++i) {
            red_indexes_values.at(s_set.at(i)) = 0;
        }
        vector<ll> result_c = computeWEFNaive(n, s, dynamic_constraints, red_indexes_values);
        ll constant = pow(2, s_set.size());
        for (int i = 0; i < result_c.size(); ++i) {
            result_c.at(i) *= constant;
        }
        result_b = add(result_b, result_c);
        return result_b;
    }
}

vector<vector<int>> get_monomials(int n) {
    vector<vector<int>> monomials;
    for (int i = 0; i < (1 << n); ++i) {
        int curRow = i;
        vector<int> b(n, 0);
        int curBinInd = 0;
        while (curRow > 0) {
            b.at(curBinInd) = curRow % 2;
            curBinInd++;
            curRow /= 2;
        }
        vector<int> monomial(n, 0);
        for (int j = 0; j < n; ++j) {
            monomial.at(j) = 1 - b.at(n - 1 - j);
        }
        monomials.push_back(monomial);
    }
    return monomials;
}


int find_the_most_right_one_pos(std::vector<int> &row_vector) {
    int pos = row_vector.size() - 1;
    while (pos > -1) {
        if (row_vector[pos]) {
            return pos;
        }
        --pos;
    }
    return pos;
}

int generate_random_bool_vector(vector<int>& odds, mt19937& gen, int odds_size) {
    int is_not_zero = 0;
    uniform_int_distribution<> dis(0, 1);
    for (int i = 0; i < odds_size; ++i) {
        int random_bit = dis(gen);
//        cout << random_bit;
        odds.push_back(random_bit);
        if (random_bit == 1) {
            is_not_zero = 1;
        }
    }
//    cout << '\n' << flush;
    return is_not_zero;
}

void apply_random_constraints(vector<pair<bool, vector<int> > >& dynamic_constraints, vector<int> frozen_bits,
                              vector<int> rows_to_add) {
    vector<int> info_indexes;
    for (int i = 0; i < frozen_bits.size(); ++i) {
        if (frozen_bits.at(i) == 0) {
            info_indexes.push_back(i);
//            cout << info_indexes.at(info_indexes.size() - 1) << ' ';
        }
    }
    vector<vector<int>> info_set_for_frozen_bit(frozen_bits.size());
    for (int i = 0; i < frozen_bits.size(); ++i) {
        if (frozen_bits.at(i) == 1) {
            for (int j = 0; j < info_indexes.size(); ++j) {
                if (i > info_indexes.at(j)) {
                    info_set_for_frozen_bit.at(i).push_back(info_indexes.at(j));
                }
            }
        }
    }
//    for (int i = 0; i < frozen_bits.size(); ++i) {
//        cout << "\n";
//        cout << i << ": ";
//        if (!info_set_for_frozen_bit.at(i).empty()) {
//            for (int j = 0; j < info_set_for_frozen_bit.at(i).size(); ++j) {
//                cout << info_set_for_frozen_bit.at(i).at(j) << ' ';
//            }
//            cout << "\n";
//        }
//    }
    random_device rd;
    mt19937 generator(rd());
//    for (int i = 0; i < rows_to_add.size(); ++i) {
//        cout << rows_to_add.at(i) << ": ";
//        for (int j = 0; j < info_set_for_frozen_bit.at(rows_to_add.at(i)).size(); ++j) {
//            cout << info_set_for_frozen_bit.at(rows_to_add.at(i)).at(j) << " ";
//        }
//        cout << '\n';
//    }

    for (int i = 0; i < rows_to_add.size(); ++i) {
        dynamic_constraints.at(rows_to_add.at(i)).first = true;
        for (int j = 0; j < info_set_for_frozen_bit.at(rows_to_add.at(i)).size(); ++j) {
            uniform_int_distribution<> dis(0, 1);
            int random_number = dis(generator);
                if (random_number == 1) {
                    dynamic_constraints.at(rows_to_add.at(i)).second.push_back(info_set_for_frozen_bit.at(rows_to_add.at(i)).at(j));
                }
        }
    }

}

void simplify_matrix(vector<pair < bool, vector<int> > >& dc, vector<int>& frozen_bits) {
    for (int i = 0; i < dc.size(); ++i) {
        if (dc.at(i).first && !dc.at(i).second.empty()) {
            vector<int> new_c;
            for (int j = 0; j < dc.at(i).second.size(); ++j) {
                if (frozen_bits.at(dc.at(i).second.at(j)) == 0) {
                    new_c.push_back(dc.at(i).second.at(j));
                } else {
                    if (!dc.at(dc.at(i).second.at(j)).second.empty()) {

                        for (int k = 0; k < dc.at(dc.at(i).second.at(j)).second.size(); ++k) {
                            int ind = -1;
                            for (int ii = 0; ii < new_c.size(); ++ii) {
                                if (new_c.at(ii) == dc.at(dc.at(i).second.at(j)).second.at(k)) {
                                    ind = ii;
                                }
                            }
                            if (ind == -1) {
                                new_c.push_back(dc.at(dc.at(i).second.at(j)).second.at(k));
                            } else {
                                new_c.erase(new_c.begin() + ind);
                            }
                        }
                    }
                }
            }

            sort(new_c.begin(), new_c.end());
            dc.at(i).second = new_c;
        }
    }
}

int mainFunc(std::shared_ptr<mc::MATLABApplication> app, const int argc, const char * argv[])
{
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    freopen("input14.txt", "r", stdin);
//    freopen("output.txt", "w", stdout);
    auto start = high_resolution_clock::now();
    int snr_number = 0;
    cin >> snr_number;
    int iteration_number;
    cin >> iteration_number; // number of iterations
    int m, k, n;
    cin >> m >> k;
    n = (1 << m);
    int v_matrix_size;
    cin >> v_matrix_size;
    vector <vector<int>> dynamic_constraints_matrix(v_matrix_size, vector<int>(n));
    for (int i = 0; i < v_matrix_size; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> dynamic_constraints_matrix.at(i).at(j);
        }
    }
    for (int snr_elem = 0; snr_elem < snr_number; ++snr_elem) {
        int last_frozen_pos = 0;
        double SNRb;
        cin >> SNRb;
        vector<int> frozen_bits(n, 0);
        for (int i = 0; i < frozen_bits.size(); ++i) {
            cin >> frozen_bits.at(i);
            if (frozen_bits.at(i) == 1) {
                last_frozen_pos = i;
            }
        }

        vector < pair < bool, vector < int > > > dynamic_constraints_with_min_bound_probability(n, {false, {}});
        double min_bound_probability = 1;
        vector<ll> wef(n + 1, 0);
        vector<int> rows_to_add;

        vector < pair < bool, vector < int > > > base_dynamic_constraints(n, {false, {}});



        for (int i = 0; i < dynamic_constraints_matrix.size(); ++i) {
            int r_pos = find_the_most_right_one_pos(dynamic_constraints_matrix.at(i));
            base_dynamic_constraints.at(r_pos).first = true;
            for (int j = 0; j < r_pos; ++j) {
                if (dynamic_constraints_matrix.at(i).at(j) == 1) {
                    base_dynamic_constraints.at(r_pos).second.push_back(j);
                }
            }
        }
        for (int i = 0; i < frozen_bits.size(); ++i) {
            if (!base_dynamic_constraints.at(i).first && frozen_bits.at(i) == 1) {
                rows_to_add.push_back(i);
            }
        }
//    for (int i = 0; i < rows_to_add.size(); ++i) {
//
//            cout << rows_to_add.at(i) << ' ';
//
//    }
//    cout << '\n';
        for (int iter = 0; iter < iteration_number; ++iter) {

            vector < pair < bool, vector < int > > > dynamic_constraints = base_dynamic_constraints;
            vector <pair<int, bool>> red_indexes;
            vector<int> red_indexes_values(n, -2);
            vector<int> info_indexes;
            for (int i = 0; i < frozen_bits.size(); ++i) {
                if (frozen_bits.at(i) == 0) {
                    if (i < last_frozen_pos) {
                        red_indexes.push_back({i, 1});
                        red_indexes_values.at(i) = -1;
                        info_indexes.push_back(i);

                    }
                }
            }

            apply_random_constraints(dynamic_constraints, frozen_bits, rows_to_add);

//            simplify_matrix(dynamic_constraints, frozen_bits);
//            cout << "\ndynamic_constraints:\n";
//            for (int i = 0; i < dynamic_constraints.size(); ++i) {
//                if (dynamic_constraints.at(i).first && !dynamic_constraints.at(i).second.empty()) {
//                    cout << i << " ";
//                    for (int j = 0; j < dynamic_constraints.at(i).second.size(); ++j) {
//                        cout << dynamic_constraints.at(i).second.at(j) << " ";
//                    }
//                    cout << '\n';
//                }
//            }

//    vector<ll> result_large1 = computeWEFNaive(n, last_frozen_pos, dynamic_constraints, red_indexes_values);
//    vector<ll> result1 = vector<ll>(result_large1.begin(), result_large1.end() - result_large1.size() + n + 1);
//    cout << "\n";
//    for (int i = 0; i < result1.size(); ++i) {
//        cout << i << ": " << result1.at(i) << '\n';
//    }
//    cout << "\n";
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<seconds>(stop - start);
//    cout << '\n';
//    cout << "Time: " << duration.count() << '\n';

            vector <vector<int>> monomials = get_monomials(m);
            vector <vector<int>> monomials_order = get_monomials_order(monomials, info_indexes);
            vector<ll> result_large = computeWEF(n, last_frozen_pos, monomials_order,
                                                 red_indexes, dynamic_constraints, red_indexes_values);

            vector<ll> result = vector<ll>(result_large.begin(), result_large.end() - result_large.size() + n + 1);

            vector<double> result_double(result.size(), 0);
//            cout << "\n";
            for (int i = 0; i < result.size(); ++i) {
                result_double.at(i) = (double) result.at(i);
//                cout << result.at(i) << ' ';
            }
//            cout << "\n";
            std::initializer_list<double> i_list(result_double.data(), result_double.data() + result_double.size());


//            auto stop1 = high_resolution_clock::now();
//            auto duration1 = duration_cast<seconds>(stop1 - start);
//            cout << '\n';
//            cout << "Time: " << duration1.count() << '\n';

            md::ArrayFactory factory;
            md::TypedArray<double> nIn = factory.createArray<double>({1, 1}, {(double) n});
            md::TypedArray<double> kIn = factory.createArray<double>({1, 1}, {(double) k});
            md::TypedArray<double> SNRbIn = factory.createArray<double>({1, 1}, {SNRb});
            md::TypedArray<double> SPIn = factory.createArray<double>({1, result_double.size()}, i_list);
            try {
                // The path to the CTF (library archive file) passed to
                // initMATLABLibrary or initMATLABLibraryAsync may be either absolute
                // or relative. If it is relative, the following will be prepended
                // to it, in turn, in order to find the CTF:
                // - the directory named by the environment variable
                // CPPSHARED_BASE_CTF_PATH, if defined
                // - the working directory
                // - the directory where the executable is located
                // - on Mac, the directory three levels above the directory
                // where the executable is located

                // If the CTF is not in one of these locations, do one of the following:
                // - copy the CTF
                // - move the CTF
                // - change the working directory ("cd") to the location of the CTF
                // - set the environment variable to the location of the CTF
                // - edit the code to change the path
                auto lib = mc::initMATLABLibrary(app, u"libCalcPoltirevBound.ctf");
                std::vector <md::Array> inputs{nIn, kIn, SNRbIn, SPIn};
                auto result1 = lib->feval(u"calcBound", 1, inputs);
                for (auto r : result1) {
                    double res = r[0];
                    if (res < min_bound_probability) {

                        min_bound_probability = res;
                        dynamic_constraints_with_min_bound_probability = dynamic_constraints;
                        wef = result;
                    }
//                    std::cout << res << std::endl;
                }
//                auto stop2 = high_resolution_clock::now();
//
//
//                auto duration2 = duration_cast<seconds>(stop2 - start);
//                cout << '\n';
//                cout << "Time: " << duration2.count() << '\n';
            } catch (const std::exception &exc) {
                std::cerr << exc.what() << std::endl;
                return -1;
            }
        }
        auto stop2 = high_resolution_clock::now();


        auto duration2 = duration_cast<seconds>(stop2 - start);
        cout << '\n';
        cout << "Time: " << duration2.count() << '\n';
        cout << "SNR: " << SNRb << "\n";
        cout << "Min Poltirev's bound spec\n";
        cout << "================================================\n";
        cout << "\ndynamic_constraints:\n";
        for (int i = 0; i < dynamic_constraints_with_min_bound_probability.size(); ++i) {
            if (dynamic_constraints_with_min_bound_probability.at(i).first
                && !dynamic_constraints_with_min_bound_probability.at(i).second.empty()) {
                cout << i << " ";
                for (int j = 0; j < dynamic_constraints_with_min_bound_probability.at(i).second.size(); ++j) {
                    cout << dynamic_constraints_with_min_bound_probability.at(i).second.at(j) << " ";
                }
                cout << '\n';
            }
        }
        cout << "\nWEF:\n";
        for (int i = 0; i < wef.size(); ++i) {
            cout << wef.at(i) << ' ';
        }
        cout << "\n";
        cout << "Min Poltirev's bound: " << min_bound_probability << '\n';
    }
	return 0;
}

// The main routine. On the Mac, the main thread runs the system code, and
// user code must be processed by a secondary thread. On other platforms, 
// the main thread runs both the system code and the user code.

// build: mbuild -O optimizer.cpp
int main(const int argc, const char * argv[])
{
	int ret = 0;
	try {
		auto matlabApplication = setup();
		ret = mc::runMain(mainFunc, std::move(matlabApplication), argc, argv);
		// Calling reset() on matlabApplication allows the user to control
		// when it is destroyed, which automatically cleans up its resources.
		// Here, the object would go out of scope and be destroyed at the end
		// of the block anyway, even if reset() were not called.
		// Whether the matlabApplication object is explicitly or implicitly
		// destroyed, initMATLABApplication() cannot be called again within
		// the same process.
		matlabApplication.reset();
	} catch(const std::exception & exc) {
		std::cerr << exc.what() << std::endl;
		return -1;
	}
	return ret;
}

