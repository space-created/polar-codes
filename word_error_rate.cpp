#include "PolarCode.h"
#include <random>

vector<double> PolarCode::get_word_error_rate(vector<double> ebno_vec, u8 l_size,
                                              size_t min_error_amount, size_t max_amount_runs) {


    vector<double> word_error(ebno_vec.size(), 0);

    double sqr_sigma = 1.0;

    std::vector<double> p0(word_length), p1(word_length);
    double sigma_sqrt_pi = std::sqrt(sqr_sigma * 3.1415f);

    normal_distribution<double> gauss_dist(0.0f, sqr_sigma);
    std::random_device rd;
    std::mt19937 generator(rd());

    for (int ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {

        size_t num_err = 0;
        size_t num_run = 0;
        while (num_err < min_error_amount && num_run < max_amount_runs) {
//            try {
            vector<u8> coded_bits;
            vector<double> bpsk(word_length);
            vector<double> received_signal(word_length, 0);
            vector<u8> info_bits = get_random_boolean_vector(info_length);
//            vector<u8> info_bits = {1,0,1,1,0,1,0};
//            vector<u8> info_bits = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0};

            vector<double> noise(word_length, 0);

            for (u16 i = 0; i < word_length; ++i) {
                noise.at(i) = (double) gauss_dist(generator);
            }

            coded_bits = encode(info_bits);
//            cout << '\n' << "Encoded vector: ";
//            for (int i = 0; i < info_bits.size(); ++i) {
//                cout << (int) info_bits.at(i) << ' ';
//            }
//            cout << '\n';
            for (u16 i = 0; i < word_length; ++i) {
                bpsk.at(i) = 2.0f * ((double) coded_bits.at(i)) - 1.0f;
            }


            double noise_amplitude = sqrt(sqr_sigma / 2);
//            double enbo_linear = pow(10.0f, ebno_vec.at(ebno_i) / 10);
//            double rate = ((double) info_length) / ((double) (word_length));
//            double signal_amplitude = noise_amplitude * sqrt(2 * rate * enbo_linear);

            double signal_amplitude = pow(10.0f, ebno_vec.at(ebno_i) / 20)
                                      * sqrt(sqr_sigma * ((double) info_length) / ((double) (word_length)));

            for (u16 i = 0; i < word_length; ++i) {
                received_signal.at(i) = signal_amplitude * bpsk.at(i) + noise_amplitude * noise.at(i);
            }

            for (u16 i = 0; i < word_length; ++i) {

                p0.at(i) =
                        exp(-(received_signal.at(i) + signal_amplitude) * (received_signal.at(i) + signal_amplitude) /
                            sqr_sigma) / sigma_sqrt_pi;
                p1.at(i) =
                        exp(-(received_signal.at(i) - signal_amplitude) * (received_signal.at(i) - signal_amplitude) /
                            sqr_sigma) / sigma_sqrt_pi;

            }

//            try {
                vector<u8> decoded_info_bits = decode(p1, p0, l_size);
//                cout << '\n' << "Decoded vector: ";
//                for (int i = 0; i < decoded_info_bits.size(); ++i) {
//                    cout << (int) decoded_info_bits.at(i) << ' ';
//                }
//                cout << '\n';
//            } catch (...) {
//                cout << "Decode " << '\n';
//            }
            num_run++;

            for (u16 i = 0; i < info_length; ++i) {
                if (info_bits.at(i) != decoded_info_bits.at(i)) {
                    num_err++;
                    break;
                }
            }
//            } catch (...) {
//                cout << "Encode " << '\n';
//            }
        }
        word_error.at(ebno_i) = (double) num_err / num_run;
    }
    return word_error;
}

