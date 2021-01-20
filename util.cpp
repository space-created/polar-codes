#include "PolarCode.h"
#include <ctime>

vector<u8> PolarCode::get_random_boolean_vector(size_t size) {
    vector<u8> random_vector(size);
    for (size_t i = 0; i < size; ++i) {
        random_vector.at(i) = rand() % 2;
    }
    return random_vector;
}