#include <iomanip>
#include <chrono>

#include "PolarCode.h"

using namespace std::chrono;

int main(int argc, char *argv[]) {

    srand(time(nullptr));
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);

    auto start = high_resolution_clock::now();
//    freopen("out.txt", "w", stdout);
//    freopen("in.txt", "r", stdin);
    vector<int> distances = {28};
//    vector<int> distances = {11, 12, 16, 21};
    for (int dist = 0; dist < 1; ++dist) {
        u8 n = 10;
        u16 info_length = 512;

//    u8 n = 10;
//    u16 info_length = 512;
//    cout << "AWGN Subcode, q=64, 1024,512,28 L = 1,..16,32 1000000 iters\n";

//    u8 n = 6;
//    u16 info_length = 30;

//    u8 n = 4;
//    u16 info_length = 6;

        u16 crc_size = 0;

        // is subcode ?
        bool is_subcode = true;
        int q = 0;

        bool is_BEC = false; // if false uses AWGN
        double epsilon = 0.5;

        vector<u16> poly(n + 1);
//            x^11 + x^2 + 1
//        poly = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1};
//            x^11 + x^5 + x^3 + x^1 + 1
//        poly = {1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1};
        u16 bch_code_distance = distances[dist];
//            x^10 + x^3 + 1
    poly = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
//    u16 bch_code_distance = 28;

//          x^6+x^1+1
//    poly = {1,0,0,0,0,1,1};
//    u16 bch_code_distance = 14;

//           x^4+x^3+1
//    poly = {1,1,0,0,1};
//    u16 bch_code_distance = 6;

        const double sigma_sqr = 1.0;

        PolarCode polar_code(n, info_length, is_BEC, epsilon, crc_size, is_subcode, poly, bch_code_distance, q,
                             sigma_sqr);

        double ebno_log_min = 2.00;
        double ebno_log_max = 2.01;
        double ebno_log_increment = 0.25;
        vector<double> ebno_vec;

        for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment) {
            ebno_vec.push_back(ebno_log);
        }

        vector<int> list_size_arr = {
//            1,
//            2,
//            4,
//            8,
//            16,
            32,
//            64,
//            256,
        };

        const size_t min_error_amount = 100;
        const size_t max_runs_amount = 10000;


        vector<vector<double> > word_error_rate(list_size_arr.size());

        word_error_rate = polar_code.get_word_error_rate(ebno_vec, list_size_arr, min_error_amount, max_runs_amount,
                                                         sigma_sqr);
        cout << distances[dist] << '\n';
        for (size_t ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {
            cout << fixed << setprecision(2) << ebno_vec.at(ebno_i) << "\t\t";
            for (size_t i = 0; i < list_size_arr.size(); ++i) {
                cout << fixed << setprecision(10) << word_error_rate.at(i).at(ebno_i) << "\t\t" << flush;
            }
            cout << "\n";
        }
    }

    auto stop = high_resolution_clock::now();


    auto duration = duration_cast<seconds>(stop - start);
    cout << '\n';
    cout << "Time: " << duration.count() << '\n';
    return 0;
}

