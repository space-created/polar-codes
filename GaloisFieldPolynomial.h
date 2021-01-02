#include <utility>

#ifndef POLARCODE_GALOISFIELDPOLYNOMIAL_H
#define POLARCODE_GALOISFIELDPOLYNOMIAL_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <cmath>


#define u8 uint8_t
#define u16 uint16_t

using namespace std;

class GaloisFieldPolynomial {
public:

    GaloisFieldPolynomial(size_t representation, size_t poly_size) {
        this->len = poly_size;
        this->representation = representation;
        fill_polynomial();
    }

    explicit GaloisFieldPolynomial(vector<u8> poly, size_t poly_size) {
        this->len = poly_size;
        this->poly = std::move(poly);
        this->representation = find_representation();
    }


    static GaloisFieldPolynomial get_remainder_by_modulo(GaloisFieldPolynomial &a, GaloisFieldPolynomial& mod_poly);

    static GaloisFieldPolynomial multiply(GaloisFieldPolynomial& a, GaloisFieldPolynomial& b);
    static GaloisFieldPolynomial power_by_modulo(size_t deg, GaloisFieldPolynomial &power_poly, GaloisFieldPolynomial& a);

    static size_t get_non_zero_poly_pos(GaloisFieldPolynomial& a);

    size_t get_representation() {
        return representation;
    }

    size_t get_size() {
        return len;
    }

    u8 get_poly(size_t i) {
        return poly.at(i);
    }

private:
    vector<u8> poly;
    size_t representation;
    size_t len;

    void fill_polynomial() {
        poly.resize(len, 0);
        size_t temp = representation;
        for (size_t i = 0; i < len; ++i) {
            poly.at(len - i - 1) = temp % 2;
            temp /= 2;
        }
    }

    size_t find_representation() {
        size_t rep = 0;
        for (size_t i = 0; i < len; ++i) {
            rep += poly.at(i) * ((size_t ) pow(2, len - i - 1));
        }
        return rep;
    }
};


#endif //POLARCODE_GALOISFIELDPOLYNOMIAL_H
