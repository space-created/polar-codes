#include "PolarCode.h"

void PolarCode::initialize_frozen_bits() {
    vector<double> channel_vec(word_length);

    for (u16 i = 0; i < word_length; ++i) {
        channel_vec.at(i) = epsilon;
    }
    for (u8 it = 0; it < m; ++it) {
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
    cout << '\n';
    for (int i = 0; i < channel_vec.size(); ++i) {
        cout << channel_vec.at(i) << ' ';
    }
    cout << '\n';
    cout << '\n';

    channel_order_descending.resize(word_length);
    size_t n_t(0);
    generate(begin(channel_order_descending), end(channel_order_descending), [&] { return n_t++; });

    sort(begin(channel_order_descending),
         end(channel_order_descending),
         [&](int i1, int i2) {
//             return channel_bit_reverse_vec.at(i1) < channel_bit_reverse_vec.at(i2);
             return channel_vec.at(bit_rev_matrix_order.at(i1)) < channel_vec.at(bit_rev_matrix_order.at(i2));
         });
    u16 effective_info_length = info_length;
    cout << "channel_order_descending: ";
    for (int i = 0; i < channel_order_descending.size(); i++) {
        cout << channel_order_descending.at(i) << ' ';
    }
    cout << '\n';
    if (is_subcode) {
        for (int i = 0; i < J.size(); ++i) {
            frozen_bits.at(J.at(i)) = 1;
        }
        int amount_to_freeze = word_length - info_length - J.size();
        int pos = channel_order_descending.size() - 1;
        while (amount_to_freeze > 0) {
            if (frozen_bits.at(channel_order_descending.at(pos)) != 1) {
                frozen_bits.at(channel_order_descending.at(pos)) = 1;
                amount_to_freeze--;
            }
            pos--;
        }
    } else {
        for (u16 i = 0; i < effective_info_length; ++i) {
            frozen_bits.at(channel_order_descending.at(i)) = 0;
        }
        for (u16 i = effective_info_length; i < word_length; ++i) {
            frozen_bits.at(channel_order_descending.at(i)) = 1;
        }
    }
    cout << '\n'<< "frozen_bits: ";
    for (int i = 0; i < frozen_bits.size(); ++i) {
        cout << (int) frozen_bits.at(i) << ' ';
    }
    cout << '\n';
    crc_matrix.resize(crc_size);
    for (u8 bit = 0; bit < crc_size; ++bit) {
        crc_matrix.at(bit).resize(info_length);
        crc_matrix.at(bit) = get_random_boolean_vector(info_length);
    }
}

void PolarCode::get_bit_rev_order() {
    for (u16 i = 0; i < word_length; ++i) {
        u16 num_to_be_reversed = i;
        bit_rev_matrix_order.at(i) = (u16) ((num_to_be_reversed & 1) << (m - 1));
        for (u8 j = (u8) (m - 1); j > 0; --j) {
            num_to_be_reversed >>= 1;
            bit_rev_matrix_order.at(i) += (num_to_be_reversed & 1) << (j - 1);
        }
    }
}