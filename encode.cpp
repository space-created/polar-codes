#include "PolarCode.h"

vector<u8> PolarCode::encode(vector<u8> info_bits) {

    vector<u8> info_bits_padded(word_length, 0);
    vector<u8> coded_bits(word_length);

    // place info bits in positions of the most strong channels
    if (is_subcode) {
        int i = 0;
        int pos = 0;
        while (i < info_length) {
            u8 info_bit = info_bits.at(i);
            if (frozen_bits.at(channel_order_descending.at(pos)) != 1) {
                info_bits_padded.at(channel_order_descending.at(pos)) = info_bit;
                i++;
            }
            pos++;
        }
    } else {
        for (size_t i = 0; i < info_length; ++i) {
            info_bits_padded.at(channel_order_descending.at(i)) = info_bits.at(i);
        }
    }



    if (is_subcode) {
        for (int i = J.size() - 1; i >= 0; --i) {
            int right_one_pos = J.at(i);
            int row_num = i;
            for (int s = 0; s < right_one_pos; ++s) {


                info_bits_padded.at(right_one_pos) =
                        (info_bits_padded.at(right_one_pos) +
                         info_bits_padded.at(s) * constraint_matrix.at(row_num).at(s)) % 2;
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

    for (size_t iteration = 0; iteration < m; ++iteration) {
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