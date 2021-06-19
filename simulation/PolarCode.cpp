#include "PolarCode.h"

void PolarCode::initialize_channel_order() {
    vector<double> channel_vec;
    if (is_BEC) {
        channel_vec.resize(word_length);
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
    } else {
        // the SIMPLIFIED GAUSSIAN APPROXIMATION was used here (http://dcn.icc.spbstu.ru/~petert/papers/rchain.pdf p.3)
        channel_vec.resize(1, 2 / sigma_sqr);
        for (int it = 0; it < m; ++it) {
            vector<double> temp = channel_vec;
            channel_vec.resize(2 * temp.size(), 0);
            for (int i = 0; i < temp.size(); ++i) {
                if (temp.at(i) > 12) {
                    channel_vec.at(2 * i) = 0.9861 * temp.at(i) - 2.3152;
                } else if (temp.at(i) > 3.5 && temp.at(i) <= 12) {
                    channel_vec.at(2 * i) = temp.at(i) * (9.005 * 0.001 * temp.at(i) + 0.7694) - 0.9507;
                } else if (temp.at(i) > 1 && temp.at(i) <= 3.5) {
                    channel_vec.at(2 * i) = temp.at(i) * (0.062883 * temp.at(i) + 0.3678) - 0.1627;
                } else {
                    channel_vec.at(2 * i) = temp.at(i) * (0.2202 * temp.at(i) + 0.06448);
                }

                channel_vec.at(2 * i + 1) = temp.at(i) + temp.at(i);
            }
        }
    }
//    cout << '\n' << channel_vec.size() << '\n';
//    for (int i = 0; i < channel_vec.size(); ++i) {
//        cout << channel_vec.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n';

    channel_order_descending.resize(word_length);
    size_t n_t(0);
    generate(begin(channel_order_descending), end(channel_order_descending), [&] { return n_t++; });

    sort(begin(channel_order_descending),
         end(channel_order_descending),
         [&](int i1, int i2) {
             if (is_BEC) {
                 return channel_vec.at(bit_rev_matrix_order.at(i1)) < channel_vec.at(bit_rev_matrix_order.at(i2));
             } else {
                 return channel_vec.at(i1) > channel_vec.at(i2);
             }
         });
//    cout << channel_order_descending.size() << '\n';
//    cout << "channel_order_descending: ";
//    for (int i = 0; i < channel_order_descending.size(); i++) {
//        cout << channel_order_descending.at(i) << ' ';
//    }
//    cout << '\n';
}

void PolarCode::initialize_frozen_bits() {
    u16 effective_info_length = info_length;
    if (is_subcode) {
        for (int i = 0; i < J.size(); ++i) {
            frozen_bits.at(J.at(i)) = 1;
        }
        int amount_to_freeze = word_length - info_length - J.size();
        int pos = channel_order_descending.size() - 1;
        while (amount_to_freeze > 0) {
            if (frozen_bits.at(channel_order_descending.at(pos)) != 1) {
                frozen_bits.at(channel_order_descending.at(pos)) = 1;
                static_frozen_channels.push_back(channel_order_descending.at(pos));
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
//    cout << '\n';
//    cout << '\n'<< "static_frozen_channels: ";
//    for (int i = 0; i < static_frozen_channels.size(); ++i) {
//        cout << (int) static_frozen_channels.at(i) << ' ';
//    }
//    cout << '\n';
//    cout << '\n' << flush;
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