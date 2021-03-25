#include "PolarCode.h"
#include <random>
#include <chrono>

using namespace std::chrono;

vector<vector<double> > PolarCode::get_word_error_rate(const vector<double>& ebno_vec, vector<int> list_size_arr,
                                                       size_t min_error_amount, size_t max_amount_runs, double sqr_sigma) {
    auto start = high_resolution_clock::now();

    vector<vector<int> > num_run(list_size_arr.size(), vector<int>(ebno_vec.size(), 0));
    vector<vector<double> > num_error(list_size_arr.size(), vector<double>(ebno_vec.size(), 0));
    vector<vector<double> > word_error_rate(list_size_arr.size(), vector<double>(ebno_vec.size(), 0));


    vector<u8> coded_bits;
    vector<double> bpsk(word_length);
    vector<double> received_signal(word_length, 0);


    std::vector<double> p0(word_length), p1(word_length);
    double sigma_sqrt_pi = sqrt(sqr_sigma * 3.1415f);

    normal_distribution<double> gauss_dist(0.0f, sqr_sigma);

    std::random_device rd;
    std::mt19937 generator(rd());

    vector<double> noise(word_length, 0);
    vector<u8> info_bits(info_length, 0);

    for (int run = 0; run < max_amount_runs; ++run) {
        if ((run % (max_amount_runs / 100)) == 0) {
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(t2 - start).count();
            cout << "Iteration: " << run << "; time: " << duration / 1000 / 1000 << " seconds; complete: "
                 << (100 * run) / max_amount_runs << "." << '\n' << flush;
        }

        if ((run % 100) == 0) {
            info_bits = get_random_boolean_vector(info_length);
        }

        for (u16 i = 0; i < word_length; ++i) {
            noise.at(i) = (double) gauss_dist(generator);
        }

        coded_bits = encode(info_bits);
        for (u16 i = 0; i < word_length; ++i) {
            bpsk.at(i) = 2.0f * ((double) coded_bits.at(i)) - 1.0f;
        }

        for (int l_cur_size_index = 0; l_cur_size_index < list_size_arr.size(); ++l_cur_size_index) {
            vector<bool> prev_decoded;
            prev_decoded.resize(ebno_vec.size(), false);

            for (int ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {

                if (num_error.at(l_cur_size_index).at(ebno_i) > min_error_amount) {
                    continue;
                }

                num_run.at(l_cur_size_index).at(ebno_i)++;

                bool run_sim = true;

                for (unsigned i_ebno2 = 0; i_ebno2 < ebno_i; ++i_ebno2) {
                    if (prev_decoded.at(i_ebno2)) {
                        run_sim = false;
                    }
                }

                if (!run_sim) {
                    continue;
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
                            exp(-(received_signal.at(i) + signal_amplitude) *
                                (received_signal.at(i) + signal_amplitude) /
                                sqr_sigma) / sigma_sqrt_pi;
                    p1.at(i) =
                            exp(-(received_signal.at(i) - signal_amplitude) *
                                (received_signal.at(i) - signal_amplitude) /
                                sqr_sigma) / sigma_sqrt_pi;

                }


                vector<u8> decoded_info_bits = decode(p1, p0, list_size_arr.at(l_cur_size_index));

                bool err = false;
                for (uint16_t i = 0; i < info_length; ++i) {
                    if (info_bits.at(i) != decoded_info_bits.at(i)) {
                        err = true;
                        break;
                    }
                }

                if (err) {
                    num_error.at(l_cur_size_index).at(ebno_i)++;
                } else {
                    prev_decoded.at(ebno_i) = true;
                }
            }
        }
    }

    for (int l_index = 0; l_index < list_size_arr.size(); ++l_index) {
        for (int i_ebno = 0; i_ebno < ebno_vec.size(); ++i_ebno) {
            word_error_rate.at(l_index).at(i_ebno) = num_error.at(l_index).at(i_ebno) / num_run.at(l_index).at(i_ebno);
        }
    }
    return word_error_rate;
}

