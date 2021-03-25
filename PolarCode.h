#include <utility>

#ifndef POLARC_POLARCODE_H
#define POLARC_POLARCODE_H

#define u8 uint8_t
#define u16 uint16_t

#include <cstdint>
#include <vector>
#include <cmath>
#include <stack>
#include <limits>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <unordered_map>

using namespace std;

class PolarCode {

public:
    PolarCode(u8 num_layers,
              u16 info_length,
              bool is_BEC,
              double epsilon,
              u16 crc_size,
              bool is_subcode,
              vector<u8> poly,
              u16 bch_distance,
              int q,
              double sigma_sqr)
            : m(num_layers),
              info_length(info_length),
              is_BEC(is_BEC),
              epsilon(epsilon),
              crc_size(crc_size),
              is_subcode(is_subcode),
              poly(std::move(poly)),
              bch_distance(bch_distance),
              q(q),
              sigma_sqr(sigma_sqr) {

        word_length = (u16) (1 << m);
        frozen_bits.resize(word_length, 0);
        bit_rev_matrix_order.resize(word_length);
        get_bit_rev_order();
        initialize_channel_order();
        if (is_subcode) {
            build_constraint_matrix();
        } else {
            initialize_frozen_bits();
        }

        list_size = 1; // default
        decoded_info_bits.resize(info_length);
    }

    vector<u8> encode(vector<u8> info_bits);

    vector<u8> decode(vector<double> &p1, vector<double> &p0, u16 ls);

    vector<vector<double> >
    get_word_error_rate(const vector<double>& ebno_vec, vector<int> list_size_arr, size_t min_error_amount,
                        size_t max_amount_runs, double sqr_sigma);

private:

    u8 m;
    u16 info_length;
    u16 word_length;
    u16 crc_size;

    bool is_BEC;
    double epsilon;
    double sigma_sqr;

    vector<int> frozen_bits;
    vector<u16> channel_order_descending;
    vector<int> static_frozen_channels;
    vector<vector<u8>> crc_matrix;
    vector<u16> bit_rev_matrix_order;

    void initialize_channel_order();

    void initialize_frozen_bits();

    void get_bit_rev_order();

    u16 list_size;

    // subcode:
    bool is_subcode;
    int q;
    vector<u8> poly;
    u16 bch_distance;
    vector<int> T_arr;
    vector<int> J;
    vector<vector<u8> > constraint_matrix;

    void build_constraint_matrix();

    void insert_least_reliable_rows_for_freezing_bits();

    vector<vector<u8> > build_bch_check_matrix();

    vector<vector<u8> > kronecker_product();

    vector<vector<u8> > get_a_matrix(vector<vector<u8> > &f_matrix);

    vector<vector<u8> > transpose_matrix(vector<vector<u8> > &a_matrix);

    vector<vector<u8> > get_matrix_product(vector<vector<u8> > &a, vector<vector<u8> > &b);

    void remove_zero_rows(std::vector<std::vector<u8> > &matrix);


    int find_the_most_right_one_pos(std::vector<u8> &row_vector);

    size_t right_col_max(const std::vector<std::vector<u8> > &matrix, size_t row, size_t col);

    int count_right_zero_columns(std::vector<std::vector<u8> > &matrix, int row_pos, int pos);

    void right_triangulation(std::vector<std::vector<u8> > &matrix);

    void matrix_reduction(vector<vector<u8>> &matrix);

    void sort_by_right_one(vector<vector<u8>> &matrix);

    stack<u16> inactive_path_indices;
    vector<u16> active_path;
    vector<vector<double *>> array_pointer_P;
    vector<vector<u8 *>> array_pointer_C;
    vector<u8 *> array_pointer_info;
    vector<vector<u16>> path_index_to_array_index;
    vector<stack<u16>> inactive_array_indices;
    vector<vector<u16>> array_reference_count;

    vector<u8> decoded_info_bits;

    void initialize_data_structures();

    u16 assign_initial_path();

    u16 clone_path(u16 l);

    void kill_path(u16 l);

    double *get_array_pointer_P(u16 lambda, u16 l);

    u8 *get_array_pointer_C(u16 lambda, u16 l);

    void recursively_calc_P(u16 lambda, u16 phi);

    void recursively_update_C(u16 lambda, u16 phi);

    void continue_paths_frozen_bit(u16 phi);

    void continue_paths_unfrozen_bit(u16 phi);

    u16 find_most_probable_path(bool check_crc);

    bool crc_check(const u8 *info_bits_padded);

    static vector<u8> get_random_boolean_vector(size_t size);
};


#endif //POLARC_POLARCODE_H
