#include "PolarCode.h"

vector<u8> PolarCode::encode(vector<u8> info_bits) {

    vector<u8> info_bits_padded(word_length, 0);
    vector<u8> coded_bits(word_length);

    for (size_t i = 0; i < info_length; ++i) {
        info_bits_padded.at(channel_order_descending.at(i)) = info_bits.at(i);
    }
//    for (int i = 0; i < channel_order_descending.size(); ++i) {
//        cout << channel_order_descending.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';

//    for (int i = 0; i < frozen_bits.size(); ++i) {
//        cout << frozen_bits.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';
//    for (int i = 0; i < info_bits.size(); ++i) {
//        cout << (int) info_bits.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';
//    for (int i = 0; i < info_bits_padded.size(); ++i) {
//        cout << (int) info_bits_padded.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';

    if (is_subcode) {
//        vector<pair<int, int> > j_i(word_length - info_length);
//        for (size_t i = 0; i < j_i.size(); ++i) {
//            int the_most_right_one_pos = find_the_most_right_one_pos(constraint_matrix.at(i));
////            j_i.at(i) = {the_most_right_one_pos, i};
//            frozen_bits_num_map.at(the_most_right_one_pos) = i;
//        }
//        sort(begin(j_i),
//             end(j_i),
//             [&](pair<int, int> p1, pair<int, int> p2) {
//                 return p1.first < p2.first;
//             });
//        for (size_t i = 0; i < j_i.size(); ++i) {
//            int right_one_pos = j_i.at(i).first;
//            int row_num = j_i.at(i).second;
//            for (int s = 0; s < right_one_pos; ++s) {
//                info_bits_padded.at(right_one_pos) =
//                        (info_bits_padded.at(right_one_pos) +
//                         info_bits_padded.at(s) * constraint_matrix.at(row_num).at(s)) % 2;
//            }
//        }
//        for (int i = 0; i < frozen_bits_num_map.size(); ++i) {
//            if (frozen_bits_num_map.at(i) != -1) {
//                int right_one_pos = i;
//                int row_num = frozen_bits_num_map.at(i);
//                for (int s = 0; s < right_one_pos; ++s) {
//                    info_bits_padded.at(right_one_pos) =
//                            (info_bits_padded.at(right_one_pos) +
//                             info_bits_padded.at(s) * constraint_matrix.at(row_num).at(s)) % 2;
//                }
//            }
//        }
        for (int i = J.size() - 1; i >= 0; --i) {
            int right_one_pos = J.at(i);
            int row_num = i;
//            cout << "\nu[" + to_string(J.at(i)) + "] = 0";
            for (int s = 0; s < right_one_pos; ++s) {
//                if (constraint_matrix.at(row_num).at(s) != 0 && info_bits_padded.at(s) != 0) {
//
//                    cout << " + u[" + to_string(s) + "]";
//                }

                info_bits_padded.at(right_one_pos) =
                        (info_bits_padded.at(right_one_pos) +
                         info_bits_padded.at(s) * constraint_matrix.at(row_num).at(s)) % 2;
            }

        }
    }
//    cout << '\n';
//    cout << '\n';
//    for (int i = 0; i < info_bits_padded.size(); ++i) {
//        cout << (int) info_bits_padded.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';
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