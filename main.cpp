#include <iomanip>
#include <chrono>

#include "PolarCode.h"

using namespace std::chrono;

int main(int argc, char *argv[]) {

    srand(time(nullptr));
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    auto start = high_resolution_clock::now();
//    freopen("out3.txt", "w", stdout);
    freopen("matrix1.txt", "r", stdin);

/*
 * Considered modulation channel
 *
 * if *is_BEC* is true: the binary erasure channel is used
 * otherwise: the additive white gaussian noise channel is used
 */
    bool is_BEC = false;
    double epsilon_BEC = 0.5;

/*
 * Apply constraints (polar subcodes)
 */
    bool apply_constraints = true;
/*
 * Number of dynamic constraints
 */
    vector<int> vec_q = {-1};

/*
 * Specs
 */

    u8 m = 11;
    u16 info_length = 1536; // bch_code_distance = 12, 24
    vector<u16> poly = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1};
//                      x^11 + x^2 + 1

//    u8 m = 10;
//    u16 info_length = 512; // bch_code_distance = 28
//    vector<u16> poly = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
//                      x^10 + x^3 + 1

//    u8 m = 6;
//    u16 info_length = 36; // bch_code_distance = 14
//    vector<u16> poly = {1, 0, 0, 0, 0, 1, 1};
//                      x^6 + x^1 + 1


//    u8 m = 4;
//    u16 info_length = 6; // bch_code_distance = 6
//    vector<u16> poly = {1, 1, 0, 0, 1};
//                      x^4 + x^3 + 1


/*
 * Distances considered
 */
//    11, 12, 16, 21, 24
    vector<int> distances = {12};

/*
 * Signal to noise ratio considered
 */
    double ebno_log_min = 1.50;
    double ebno_log_max = 1.51;
    double ebno_log_increment = 0.25;

    vector<double> ebno_vec;

    for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment) {
        ebno_vec.push_back(ebno_log);
    }
/*
 * Simulation runs
 */
    const size_t min_error_amount = 100;
    const size_t max_runs_amount = 10000;
/*
 * List sizes considered
 */
    vector<int> list_size_arr = {
//            1,
//            2,
//            4,
//            8,
            16,
            32,
//            64,
//            256,
    };

/*
 * cyclic redundancy check bits number
 */
    u16 crc_size = 0;


    for (int q_num = 0; q_num < vec_q.size(); ++q_num) {
        u16 word_length = (1 << m);
        cout << "Modulation channel: " << (is_BEC ? "BEC" : "AWGN") << "\n";
        cout << (apply_constraints ? "Polar subcode (" : "Arikan Polar code (");
        cout << word_length << ", " << info_length;
        if (apply_constraints) {
            cout << "), q = " << vec_q.at(q_num) << "\n";
        } else {
            cout << ")\n";
        }

        /*
         * Main loop
         */
        for (int dist = 0; dist < distances.size(); ++dist) {
            cout << "Distance: " << distances.at(dist) << "\n";
            u16 bch_code_distance = distances.at(dist);
            cout << "SNR \t\tNoiseStdDev \t\t";
            for (int list_size_num = 0; list_size_num < list_size_arr.size(); ++list_size_num) {
                cout << list_size_arr.at(list_size_num) << "         \t\t";
            }
            cout << "\n";

            for (int ebno_num = 0; ebno_num < ebno_vec.size(); ++ebno_num) {
                cout << fixed << setprecision(2) << ebno_vec.at(ebno_num) << "\t\t";
                double noiseStdDev = sqrt(
                        0.5 * pow(10, -ebno_vec.at(ebno_num) / 10) / ((double) info_length / (double) (word_length)));
                cout << fixed << setprecision(11) << noiseStdDev << "\t\t" << flush;

                PolarCode polar_code(m, info_length, is_BEC, epsilon_BEC, crc_size, apply_constraints, poly,
                                     bch_code_distance, vec_q.at(q_num),
                                     noiseStdDev * noiseStdDev);

                for (int l_cur_size_index = 0; l_cur_size_index < list_size_arr.size(); ++l_cur_size_index) {
                    double word_error_rate = polar_code.get_word_error_rate(
                            noiseStdDev,
                            list_size_arr.at(l_cur_size_index),
                            min_error_amount,
                            max_runs_amount);
                    cout << fixed << setprecision(10) << word_error_rate << "\t\t" << flush;
                }
                cout << "\n";
            }
        }
    }
    auto stop = high_resolution_clock::now();


    auto duration = duration_cast<seconds>(stop - start);
    cout << '\n';
    cout << "Time: " << duration.count() << '\n';
    return 0;
}

