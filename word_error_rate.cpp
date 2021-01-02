#include "PolarCode.h"
#include <random>

vector<double> PolarCode::get_word_error_rate(vector<double> ebno_vec, u8 l_size,
                                              size_t min_error_amount, size_t max_amount_runs) {


    vector<double> word_error(ebno_vec.size(), 0);

    double N_0 = 1.0;

    std::vector<double> p0(word_length), p1(word_length);
    double sigma_sqrt_pi = std::sqrt(N_0 * 3.1415f);

    normal_distribution<double> gauss_dist(0.0f, N_0);
    std::random_device rd;
    std::mt19937 generator(rd());

    for (int ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {
        size_t num_err = 0;
        size_t num_run = 0;
        while (num_err < min_error_amount && num_run < max_amount_runs) {
            vector<u8> coded_bits;
            vector<double> bpsk(word_length);
            vector<double> received_signal(word_length, 0);
            vector<u8> info_bits(info_length, 0);

            vector<double> noise(word_length, 0);

            for (u16 i = 0; i < info_length; ++i) {
                info_bits.at(i) = (u8) (get_random_bool() % 2);
            }

            for (u16 i = 0; i < word_length; ++i) {
                noise.at(i) = (double) gauss_dist(generator);
            }

            coded_bits = encode(info_bits);

            for (u16 i = 0; i < word_length; ++i) {
                bpsk.at(i) = 2.0f * ((double) coded_bits.at(i)) - 1.0f;
            }


            double snr_sqrt_linear = pow(10.0f, ebno_vec.at(ebno_i) / 20)
                                     * sqrt(((double) info_length) / ((double) (word_length)));
            for (u16 i = 0; i < word_length; ++i) {
                received_signal.at(i) = snr_sqrt_linear * bpsk.at(i) + sqrt(N_0 / 2) * noise.at(i);
            }
            for (u16 i = 0; i < word_length; ++i) {

                p0.at(i) =
                        exp(-(received_signal.at(i) + snr_sqrt_linear) * (received_signal.at(i) + snr_sqrt_linear) /
                            N_0) / sigma_sqrt_pi;
                p1.at(i) =
                        exp(-(received_signal.at(i) - snr_sqrt_linear) * (received_signal.at(i) - snr_sqrt_linear) /
                            N_0) / sigma_sqrt_pi;

            }
            vector<u8> decoded_info_bits = decode(p1, p0, l_size);
            num_run++;

            for (u16 i = 0; i < info_length; ++i) {
                if (info_bits.at(i) != decoded_info_bits.at(i)) {
                    num_err++;
                    break;
                }
            }
        }
        word_error.at(ebno_i) = (double) num_err / num_run;
    }

    return word_error;
}

