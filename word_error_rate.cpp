#include "PolarCode.h"
#include <random>
#include <chrono>

#define PI 3.1415926535f

using namespace std::chrono;

double PolarCode::get_word_error_rate(double noiseStdDev, int ls,
                                      size_t min_error_amount, size_t max_amount_runs) {
    list_size = ls;
    auto start = high_resolution_clock::now();

    double num_error = 0;
    int num_run = 0;

    vector<u8> coded_bits;
    vector<double> bpsk(word_length);
    vector<double> received_signal(word_length, 0);


    std::vector<double> p0(word_length), p1(word_length);

    std::random_device rd;
    std::mt19937 generator(rd());

    vector<double> noise(word_length, 0);
    vector<u8> info_bits(info_length, 0);

    for (int run = 0; run < max_amount_runs; ++run) {
//        if ((run % (max_amount_runs / 100)) == 0) {
//            high_resolution_clock::time_point t2 = high_resolution_clock::now();
//            auto duration = duration_cast<microseconds>(t2 - start).count();
//            cout << "\nIteration: " << run << "; time: " << duration / 1000 / 1000 << " seconds; complete: "
//                 << (100 * run) / max_amount_runs << "." << '\n' << flush;
//        }

        info_bits = get_random_boolean_vector(info_length);

        coded_bits = encode(info_bits);
        for (u16 i = 0; i < word_length; ++i) {
            bpsk.at(i) = 2.0f * ((double) coded_bits.at(i)) - 1.0f;
        }
        double norm_constant = 1 / noiseStdDev * sqrt(2 * PI);

        normal_distribution<double> gauss_dist(0.0f, noiseStdDev);

        for (u16 i = 0; i < word_length; ++i) {
            noise.at(i) = (double) gauss_dist(generator);
        }

        if (num_error > min_error_amount) {
            break;
        }

        num_run++;

        for (u16 i = 0; i < word_length; ++i) {
            received_signal.at(i) = bpsk.at(i) + noise.at(i);
        }

        for (u16 i = 0; i < word_length; ++i) {
            p0.at(i) = norm_constant * exp(
                    -(received_signal.at(i) + 1.0) * (received_signal.at(i) + 1.0)
                    / (2 * noiseStdDev * noiseStdDev)
            );
            p1.at(i) = norm_constant * exp(
                    -(received_signal.at(i) - 1.0) * (received_signal.at(i) - 1.0)
                    / (2 * noiseStdDev * noiseStdDev)
            );

        }

        decode(p1, p0);

        for (uint16_t i = 0; i < info_length; ++i) {
            if (info_bits.at(i) != decoded_info_bits.at(i)) {
                num_error++;
                break;
            }
        }
    }
    return num_error / num_run;
}

