#include <iomanip>
#include <chrono>

#include "PolarCode.h"

using namespace std::chrono;

int main(int argc, char *argv[]) {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    freopen("results\\output.txt", "w", stdout);

    u8 n = 10;
    u16 info_length = (1 << (n - 1));
    u16 crc_size = 0;

    double epsilon = 0.5;

    // is subcode ?
    bool is_subcode = false;
    vector<u8> poly(n + 1, 0);
    poly = {1,0,1,0,1,1,0,1,0,1,1};
    u16 bch_info_length = 7;
    u16 bch_distance = 28;
    // is subcode ?


    PolarCode polar_code(n, info_length, epsilon, crc_size, is_subcode, poly, bch_info_length, bch_distance);

    double ebno_log_min = 1.00;
    double ebno_log_max = 3.01;
    double ebno_log_increment = 0.25;
    vector<double> ebno_vec;

    for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment) {
        ebno_vec.push_back(ebno_log);
    }

    const u8 list_size = 4;
    const size_t min_error_amount = 100;
    const size_t max_runs_amount = 100000;

    auto start = high_resolution_clock::now();

    vector<double> word_error_rate = polar_code.get_word_error_rate(ebno_vec, list_size, min_error_amount,
                                                                    max_runs_amount);

    for (size_t ebno_i = 0; ebno_i < ebno_vec.size(); ++ebno_i) {
        cout << fixed << setprecision(2) << ebno_vec.at(ebno_i) << "\t \t";
        cout << fixed << setprecision(10) << word_error_rate.at(ebno_i) << "\n";
    }

    auto stop = high_resolution_clock::now();


    auto duration = duration_cast<seconds>(stop - start);
    cout << '\n';
    cout << "Time: " << duration.count() << '\n';
    return 0;
}
//true 64
//1.00	 	0.2217294900
//1.25	 	0.1851851852
//1.50	 	0.1326259947
//1.75	 	0.1020408163
//2.00	 	0.0917431193
//2.25	 	0.0488997555
//2.50	 	0.0443852641
//2.75	 	0.0284981476
//3.00	 	0.0173010381
//false
//1.00	 	0.2512562814
//1.25	 	0.2044989775
//1.50	 	0.1533742331
//1.75	 	0.1075268817
//2.00	 	0.0728862974
//2.25	 	0.0572737686
//2.50	 	0.0433651344
//2.75	 	0.0282485876
//3.00	 	0.0156298843