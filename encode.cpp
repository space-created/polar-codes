#include "PolarCode.h"

vector<u8> PolarCode::encode(vector<u8> info_bits) {

    vector<u8> info_bits_padded(word_length, 0);
    vector<u8> coded_bits(word_length);

    for (size_t i = 0; i < info_length; ++i) {
        info_bits_padded.at(channel_order_descending.at(i)) = info_bits.at(i);
    }
    if (is_subcode) {
        frozen_bits_num_map.resize(frozen_bits.size(), -1);
        size_t i1 = 0;
        for (size_t j_i = 0; j_i < frozen_bits.size(); ++j_i) {
            if (frozen_bits.at(j_i)) {
                for (size_t s = 0; s < j_i; ++s) {
                    info_bits_padded.at(j_i) =
                            (info_bits_padded.at(j_i) + info_bits_padded.at(s) * constraint_matrix.at(i1).at(s)) % 2;
                }
                frozen_bits_num_map.at(j_i) = i1;
                i1++;
            }
        }
    }

    for (size_t i = info_length; i < info_length + crc_size; ++i) {
        u8 crc_bit = 0;
        for (size_t j = 0; j < info_length; ++j) {
            crc_bit = (u8) ((crc_bit + crc_matrix.at(i - info_length).at(j) * info_bits.at(j)) % 2);
        }
        info_bits_padded.at(channel_order_descending.at(i)) = crc_bit;
    }

    for (size_t iteration = 0; iteration < n; ++iteration) {
        auto increment = (u16) (1 << iteration);
        for (size_t j = 0; j < increment; j += 1) {
            for (size_t i = 0; i < word_length; i += 2 * increment) {
                info_bits_padded.at(i + j) = (u8) (
                        (info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment)) % 2);
            }
        }
    }

    for (size_t i = 0; i < word_length; ++i) {
        coded_bits.at(i) = info_bits_padded.at(bit_rev_matrix_order.at(i));
    }

    return coded_bits;
}