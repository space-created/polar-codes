#include <iomanip>
#include <chrono>

#include "PolarCode.h"

using namespace std::chrono;

int main(int argc, char *argv[]) {

    srand(time(nullptr));
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
//    freopen("results/output.txt", "w", stdout);

//    u8 n = 10;
//    u16 info_length = 512;

    u8 n = 6;
    u16 info_length = 36;

//    u8 n = 4; // 16
//    u16 info_length = 6; //

    u16 crc_size = 0;

    double epsilon = 0.5;

    // is subcode ?
    bool is_subcode = true;

    vector<u8> poly(n + 1);
//            x^10 + x^3 + 1
//    poly = {1,0,0,0,0,0,0,1,0,0,1};
//    u16 bch_code_distance = 28;

//          x^6+x^1+1
    poly = {1,0,0,0,0,1,1};
    u16 bch_code_distance = 12;

//           x^4+x^3+1
//    poly = {1,1,0,0,1};
//    u16 bch_code_distance = 6;
//     is subcode ?


    PolarCode polar_code(n, info_length, epsilon, crc_size, is_subcode, poly, bch_code_distance);

    double ebno_log_min = 2.50;
    double ebno_log_max = 2.51;
//    double ebno_log_max = 1.01;
    double ebno_log_increment = 0.25;
    vector<double> ebno_vec;

    for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment) {
        ebno_vec.push_back(ebno_log);
    }

    vector<u8> list_size_arr = {
//            1,
//            2,
//            4,
//            8,
            16,
//            32
    };

    auto start = high_resolution_clock::now();
    vector<vector<double> > word_error_rate(list_size_arr.size());

    for (u8 i = 0; i < list_size_arr.size(); ++i) {
        u8 list_size = list_size_arr[i];
        const size_t min_error_amount = 100;
        const size_t max_runs_amount = 100000;
//        const size_t max_runs_amount = 10000;


        word_error_rate[i] = polar_code.get_word_error_rate(ebno_vec, list_size, min_error_amount,
                                                                        max_runs_amount);
    }


    for (size_t ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {
        cout << fixed << setprecision(2) << ebno_vec.at(ebno_i) << "\t\t";
        for (size_t i = 0; i < list_size_arr.size(); ++i) {
            cout << fixed << setprecision(10) << word_error_rate.at(i).at(ebno_i) << "\t\t";
        }
        cout << "\n";
    }
    auto stop = high_resolution_clock::now();


    auto duration = duration_cast<seconds>(stop - start);
    cout << '\n';
    cout << "Time: " << duration.count() << '\n';
    return 0;
}

