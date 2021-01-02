#include "PolarCode.h"
#include <random>

u8 PolarCode::get_random_bool() {
    uniform_int_distribution<> bool_dist(0, 1);
    std::random_device rd;
    std::mt19937 mt(rd());
    return bool_dist(mt);
}