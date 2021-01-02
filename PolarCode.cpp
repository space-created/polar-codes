#include "PolarCode.h"

void PolarCode::initialize_frozen_bits() {
    vector<double> channel_vec(word_length);

    for (u16 i = 0; i < word_length; ++i) {
        channel_vec.at(i) = epsilon;
    }
    for (u8 it = 0; it < n; ++it) {
        u16 increment = 1 << it;
        for (u16 j = 0; j < increment; j += 1) {
            for (u16 i = 0; i < word_length; i += 2 * increment) {
                double z1 = channel_vec.at(i + j);
                double z2 = channel_vec.at(i + j + increment);
                channel_vec.at(i + j) = z1 + z2 - z1 * z2;
                channel_vec.at(i + j + increment) = z1 * z2;
            }
        }
    }

    channel_order_descending.resize(word_length);
    size_t n_t(0);
    generate(begin(channel_order_descending), end(channel_order_descending), [&] { return n_t++; });

    sort(begin(channel_order_descending),
         end(channel_order_descending),
         [&](int i1, int i2) {
             return channel_vec.at(bit_rev_matrix_order.at(i1)) < channel_vec.at(bit_rev_matrix_order.at(i2));
         });
    u16 effective_info_length = info_length;

    for (u16 i = 0; i < effective_info_length; ++i) {
        frozen_bits.at(channel_order_descending.at(i)) = 0;
    }
    for (u16 i = effective_info_length; i < word_length; ++i) {
        frozen_bits.at(channel_order_descending.at(i)) = 1;
    }

    crc_matrix.resize(crc_size);
    for (u8 bit = 0; bit < crc_size; ++bit) {
        crc_matrix.at(bit).resize(info_length);
        for (u16 info_bit = 0; info_bit < info_length; ++info_bit) {
            crc_matrix.at(bit).at(info_bit) = (u8) (get_random_bool() % 2);
        }
    }
    cout << endl;

    for (size_t i = 0; i < frozen_bits.size(); ++i) {
        cout << frozen_bits.at(i) << ' ';
    }
    cout << endl;

}

void PolarCode::get_bit_rev_order() {
    for (u16 i = 0; i < word_length; ++i) {
        u16 num_to_be_reversed = i;
        bit_rev_matrix_order.at(i) = (u16) ((num_to_be_reversed & 1) << (n - 1));
        for (u8 j = (u8) (n - 1); j > 0; --j) {
            num_to_be_reversed >>= 1;
            bit_rev_matrix_order.at(i) += (num_to_be_reversed & 1) << (j - 1);
        }
    }
}

