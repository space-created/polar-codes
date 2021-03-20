#include <stdexcept>

#include "GaloisFieldPolynomial.h"


GaloisFieldPolynomial GaloisFieldPolynomial::power_by_modulo(size_t deg, GaloisFieldPolynomial &power_poly, GaloisFieldPolynomial &mod_poly) {
    vector<u8> power_poly_by_modulo(power_poly.get_size(), 0);
    if (deg == 0) {
        return GaloisFieldPolynomial(1, power_poly.get_size());
    }
    auto temp_power_poly = power_poly;
    for (size_t i = 1; i < deg; ++i) {
        temp_power_poly = multiply(temp_power_poly, power_poly);
    }

    return get_remainder_by_modulo(temp_power_poly, mod_poly);
}

GaloisFieldPolynomial GaloisFieldPolynomial::multiply(GaloisFieldPolynomial &a, GaloisFieldPolynomial &b) {
    size_t a_size = a.get_size();
    size_t b_size = b.get_size();
    size_t prod_size = a_size + b_size - 1;

    vector<u8> prod_poly(prod_size, 0);

    for (int i = 0; i < a_size; i++) {
        for (int j = 0; j < b_size; j++) {
            prod_poly.at(i + j) ^= a.get_poly(i) * b.get_poly(j);
        }
    }
    return GaloisFieldPolynomial(prod_poly, prod_poly.size());
}

GaloisFieldPolynomial
GaloisFieldPolynomial::get_remainder_by_modulo(GaloisFieldPolynomial &a, GaloisFieldPolynomial &mod_poly) {
    size_t non_zero_pos_a = get_non_zero_poly_pos(a);
    size_t non_zero_pos_mod_poly = get_non_zero_poly_pos(mod_poly);

    if (non_zero_pos_mod_poly == -1) {
        throw std::overflow_error("Divide by zero exception");
    }

    if (non_zero_pos_a < non_zero_pos_mod_poly) {
        return a;
    }

    vector<u8> a_poly(non_zero_pos_a + 1);
    vector<u8> m_poly(non_zero_pos_mod_poly + 1);
    for (size_t i = 0; i < a_poly.size(); ++i) {
        a_poly.at(i) = a.get_poly(a.get_size() - 1 - non_zero_pos_a + i);
    }

    for (size_t i = 0; i < m_poly.size(); ++i) {
        m_poly.at(i) = mod_poly.get_poly(mod_poly.get_size() - 1 - non_zero_pos_mod_poly + i);
    }
    size_t shift = 0;
    while (shift <= a_poly.size() - m_poly.size()) {
        size_t new_shift = shift;
        for (size_t i = 0; i < m_poly.size(); ++i) {
            a_poly.at(i + shift) ^= m_poly.at(i);
        }

        for (size_t i = 0; i < m_poly.size(); ++i) {
            if (a_poly.at(i + shift) == 0) {
                new_shift++;
            } else {
                break;
            }
        }
        shift = new_shift;
    }

    vector<u8> remainder_poly(m_poly.size() - 1);
    for (size_t i = 0; i < remainder_poly.size(); ++i) {
        remainder_poly.at(remainder_poly.size() - i - 1) = a_poly.at(a_poly.size() - i - 1);
    }

    return GaloisFieldPolynomial(remainder_poly, remainder_poly.size());
}

size_t GaloisFieldPolynomial::get_non_zero_poly_pos(GaloisFieldPolynomial &a) {
    size_t i = 0;
    while (i < a.get_size() && a.get_poly(i) != 1) {
        i++;
    }
    if (i == a.get_size()) {
        return -1;
    }
    return a.get_size() - i - 1;
}





