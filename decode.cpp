#include "PolarCode.h"

vector<u8> PolarCode::decode(vector<double> &p1, vector<double> &p0) {

    initialize_data_structures();

    u16 l = assign_initial_path();

    double *p_0 = get_array_pointer_P(0, l);

    for (u16 beta = 0; beta < word_length; ++beta) {
        p_0[2 * beta] = (double) p0.at(beta);
        p_0[2 * beta + 1] = (double) p1.at(beta);
    }

    for (u16 phi = 0; phi < word_length; ++phi) {
        recursively_calc_P(m, phi);

        if (frozen_bits.at(phi) == 1) {
            continue_paths_frozen_bit(phi);
        } else {
            continue_paths_unfrozen_bit(phi);
        }
        if ((phi % 2) == 1) {
            recursively_update_C(m, phi);
        }
    }
    u16 most_probable_path = find_most_probable_path((bool) crc_size);

    u8 *c_0 = array_pointer_info.at(most_probable_path);

    if (is_subcode) {
        int i = 0;
        int pos = 0;
        while (i < info_length) {
            if (frozen_bits.at(channel_order_descending.at(pos)) != 1) {
                decoded_info_bits.at(i) = c_0[channel_order_descending.at(pos)];
                i++;
            }
            pos++;
        }
    } else {
        for (u16 beta = 0; beta < info_length; ++beta) {
            decoded_info_bits.at(beta) = c_0[channel_order_descending.at(beta)];
        }
    }

    for (u16 s = 0; s < list_size; ++s) {
        delete[] array_pointer_info.at(s);
        for (u16 lambda = 0; lambda < m + 1; ++lambda) {

            delete[] array_pointer_P.at(lambda).at(s);

            delete[] array_pointer_C.at(lambda).at(s);
        }
    }

    return decoded_info_bits;
}


void PolarCode::initialize_data_structures() {

    while (!inactive_path_indices.empty()) {
        inactive_path_indices.pop();
    }
    active_path.resize(list_size);


    array_pointer_P.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
        array_pointer_P.at(i).resize(list_size);
    }


    array_pointer_C.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
        array_pointer_C.at(i).resize(list_size);
    }

    array_pointer_info.resize(list_size);

    path_index_to_array_index.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
        path_index_to_array_index.at(i).resize(list_size);
    }

    inactive_array_indices.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
        while (!inactive_array_indices.at(i).empty()) {
            inactive_array_indices.at(i).pop();
        }
    }

    array_reference_count.resize(m + 1);
    for (int i = 0; i < m + 1; ++i) {
        array_reference_count.at(i).resize(list_size);
    }
    for (u16 s = 0; s < list_size; ++s) {
        array_pointer_info.at(s) = new u8[word_length]();
        for (u16 lambda = 0; lambda < m + 1; ++lambda) {

            array_pointer_P.at(lambda).at(s) = new double[2 * (1 << (m - lambda))]();

            array_pointer_C.at(lambda).at(s) = new u8[2 * (1 << (m - lambda))]();
            array_reference_count.at(lambda).at(s) = 0;
            inactive_array_indices.at(lambda).push(s);
        }
    }

    for (u16 l = 0; l < list_size; ++l) {
        active_path.at(l) = 0;
        inactive_path_indices.push(l);
    }
}

u16 PolarCode::assign_initial_path() {

    u16 l = inactive_path_indices.top();
    inactive_path_indices.pop();
    active_path.at(l) = 1;
    for (u16 lambda = 0; lambda < m + 1; ++lambda) {
        u16 s = inactive_array_indices.at(lambda).top();
        inactive_array_indices.at(lambda).pop();
        path_index_to_array_index.at(lambda).at(l) = s;
        array_reference_count.at(lambda).at(s) = 1;
    }
    return l;
}

u16 PolarCode::clone_path(u16 l) {
    u16 l_p = inactive_path_indices.top();
    inactive_path_indices.pop();
    active_path.at(l_p) = 1;

    for (u16 lambda = 0; lambda < m + 1; ++lambda) {
        u16 s = path_index_to_array_index.at(lambda).at(l);
        path_index_to_array_index.at(lambda).at(l_p) = s;
        array_reference_count.at(lambda).at(s)++;
    }
    return l_p;
}

void PolarCode::kill_path(u16 l) {
    active_path.at(l) = 0;
    inactive_path_indices.push(l);
    for (u16 lambda = 0; lambda < m + 1; ++lambda) {
        u16 s = path_index_to_array_index.at(lambda).at(l);
        array_reference_count.at(lambda).at(s)--;
        if (array_reference_count.at(lambda).at(s) == 0) {
            inactive_array_indices.at(lambda).push(s);
        }
    }
}

double *PolarCode::get_array_pointer_P(u16 lambda, u16 l) {
    u16 s = path_index_to_array_index.at(lambda).at(l);
    u16 s_p;
    if (array_reference_count.at(lambda).at(s) == 1) {
        s_p = s;
    } else {
        s_p = inactive_array_indices.at(lambda).top();
        inactive_array_indices.at(lambda).pop();


        copy(array_pointer_P.at(lambda).at(s), array_pointer_P.at(lambda).at(s) + (1 << (m - lambda + 1)),
             array_pointer_P.at(lambda).at(s_p));
        copy(array_pointer_C.at(lambda).at(s), array_pointer_C.at(lambda).at(s) + (1 << (m - lambda + 1)),
             array_pointer_C.at(lambda).at(s_p));

        array_reference_count.at(lambda).at(s)--;
        array_reference_count.at(lambda).at(s_p) = 1;
        path_index_to_array_index.at(lambda).at(l) = s_p;
    }
    return array_pointer_P.at(lambda).at(s_p);
}


u8 *PolarCode::get_array_pointer_C(u16 lambda, u16 l) {
    u16 s = path_index_to_array_index.at(lambda).at(l);
    u16 s_p;
    if (array_reference_count.at(lambda).at(s) == 1) {
        s_p = s;
    } else {

        s_p = inactive_array_indices.at(lambda).top();
        inactive_array_indices.at(lambda).pop();


        copy(array_pointer_P.at(lambda).at(s), array_pointer_P.at(lambda).at(s) + (1 << (m - lambda + 1)),
             array_pointer_P.at(lambda).at(s_p));

        copy(array_pointer_C.at(lambda).at(s), array_pointer_C.at(lambda).at(s) + (1 << (m - lambda + 1)),
             array_pointer_C.at(lambda).at(s_p));

        array_reference_count.at(lambda).at(s)--;
        array_reference_count.at(lambda).at(s_p) = 1;
        path_index_to_array_index.at(lambda).at(l) = s_p;

    }
    return array_pointer_C.at(lambda).at(s_p);
}

void PolarCode::recursively_calc_P(u16 lambda, u16 phi) {
    if (lambda == 0)
        return;
    u16 psi = phi >> 1;
    if ((phi % 2) == 0) {
        recursively_calc_P(lambda - 1, psi);
    }
    double sigma = 0.0f;
    for (u16 l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            continue;
        }
        double *p_lambda = get_array_pointer_P(lambda, l);
        double *p_lambda_1 = get_array_pointer_P(lambda - 1, l);

        u8 *c_lambda = get_array_pointer_C(lambda, l);
        for (u16 beta = 0; beta < (1 << (m - lambda)); ++beta) {
            if ((phi % 2) == 0) {
                p_lambda[2 * beta] = 0.5f * (p_lambda_1[2 * (2 * beta)] * p_lambda_1[2 * (2 * beta + 1)]
                                             + p_lambda_1[2 * (2 * beta) + 1] * p_lambda_1[2 * (2 * beta + 1) + 1]);
                p_lambda[2 * beta + 1] = 0.5f * (p_lambda_1[2 * (2 * beta) + 1] * p_lambda_1[2 * (2 * beta + 1)]
                                                 + p_lambda_1[2 * (2 * beta)] * p_lambda_1[2 * (2 * beta + 1) + 1]);
            } else {
                u8 u_p = c_lambda[2 * beta];
                p_lambda[2 * beta] = 0.5f * p_lambda_1[2 * (2 * beta) + (u_p % 2)] * p_lambda_1[2 * (2 * beta + 1)];
                p_lambda[2 * beta + 1] =
                        0.5f * p_lambda_1[2 * (2 * beta) + ((u_p + 1) % 2)] * p_lambda_1[2 * (2 * beta + 1) + 1];
            }
            sigma = max(sigma, p_lambda[2 * beta]);
            sigma = max(sigma, p_lambda[2 * beta + 1]);
        }
    }

    for (u16 l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            continue;
        }
        double *p_lambda = get_array_pointer_P(lambda, l);
        for (u16 beta = 0; beta < (1 << (m - lambda)); ++beta) {
            p_lambda[2 * beta] = p_lambda[2 * beta] / sigma;
            p_lambda[2 * beta + 1] = p_lambda[2 * beta + 1] / sigma;
        }
    }
}

void PolarCode::recursively_update_C(u16 lambda, u16 phi) {

    u16 psi = phi >> 1;
    for (u16 l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            continue;
        }
        u8 *c_lambda = get_array_pointer_C(lambda, l);
        u8 *c_lambda_1 = get_array_pointer_C(lambda - 1, l);
        for (u16 beta = 0; beta < (1 << (m - lambda)); ++beta) {
            c_lambda_1[2 * (2 * beta) + (psi % 2)] = (u8) ((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
            c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
        }
    }
    if ((psi % 2) == 1) {
        recursively_update_C((u16) (lambda - 1), psi);
    }
}

void PolarCode::continue_paths_frozen_bit(u16 phi) {
    for (u16 l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            continue;
        }
        u8 *c_m = get_array_pointer_C(m, l);
        if (is_subcode) {
            u8 value = 0;
            for (size_t s = 0; s < phi; ++s) {
                if (T_arr.at(phi) != -1) {
                    value = (value + array_pointer_info.at(l)[s] * constraint_matrix[T_arr.at(phi)][s]) % 2;
                }
            }
            c_m[(phi % 2)] = value;
            array_pointer_info.at(l)[phi] = value;
        } else {
            c_m[(phi % 2)] = 0;
            array_pointer_info.at(l)[phi] = 0;
        }
    }
}

void PolarCode::continue_paths_unfrozen_bit(u16 phi) {
    vector<double> probForks((unsigned long) (2 * list_size));
    vector<double> probabilities;
    vector<uint8_t> contForks((unsigned long) (2 * list_size));
    u16 i = 0;
    for (unsigned l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            probForks.at(2 * l) = -1;
            probForks.at(2 * l + 1) = -1;
        } else {

            double *p_m = get_array_pointer_P(m, l);
            probForks.at(2 * l) = p_m[0];
            probForks.at(2 * l + 1) = p_m[1];


            probabilities.push_back(probForks.at(2 * l));
            probabilities.push_back(probForks.at(2 * l + 1));

            i++;
        }
    }

    u16 rho = list_size;
    if ((2 * i) < list_size) {
        rho = (u16) 2 * i;
    }
    for (u16 l = 0; l < 2 * list_size; ++l) {
        contForks.at(l) = 0;
    }
    sort(probabilities.begin(), probabilities.end(), greater<double>());
    double threshold = probabilities.at((unsigned long) (rho - 1));
    u16 num_paths_continued = 0;

    for (u16 l = 0; l < 2 * list_size; ++l) {
        if (probForks.at(l) > threshold) {
            contForks.at(l) = 1;
            num_paths_continued++;
        }
        if (num_paths_continued == rho) {
            break;
        }
    }

    if (num_paths_continued < rho) {
        for (u16 l = 0; l < 2 * list_size; ++l) {
            if (probForks.at(l) == threshold) {
                contForks.at(l) = 1;
                num_paths_continued++;
            }
            if (num_paths_continued == rho) {
                break;
            }
        }
    }

    for (unsigned l = 0; l < list_size; ++l) {
        if (active_path.at(l) == 0) {
            continue;
        }
        if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0) {
            kill_path(l);
        }
    }

    for (unsigned l = 0; l < list_size; ++l) {
        if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0) {
            continue;
        }
        u8 *c_m = get_array_pointer_C(m, l);

        if (contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1) {

            c_m[(phi % 2)] = 0;
            u16 l_p = clone_path(l);
            c_m = get_array_pointer_C(m, l_p);
            c_m[(phi % 2)] = 1;

            copy(array_pointer_info.at(l), array_pointer_info.at(l) + phi, array_pointer_info.at(l_p));
            array_pointer_info.at(l)[phi] = 0;
            array_pointer_info.at(l_p)[phi] = 1;


        } else {
            if (contForks.at(2 * l) == 1) {
                c_m[(phi % 2)] = 0;
                array_pointer_info.at(l)[phi] = 0;


            } else {
                c_m[(phi % 2)] = 1;
                array_pointer_info.at(l)[phi] = 1;

            }
        }
    }
}

u16 PolarCode::find_most_probable_path(bool check_crc) {

    u16 l_p = 0;
    double p_p1 = 0;
    bool path_with_crc_pass = false;
    for (u16 l = 0; l < list_size; ++l) {

        if (active_path.at(l) == 0) {
            continue;
        }

        if ((check_crc) && (!crc_check(array_pointer_info.at(l)))) {
            continue;
        }

        path_with_crc_pass = true;

        u8 *c_m = get_array_pointer_C(m, l);
        double *p_m = get_array_pointer_P(m, l);
        if (p_p1 < p_m[c_m[1]]) {
            l_p = l;
            p_p1 = p_m[c_m[1]];
        }
    }
    if (path_with_crc_pass) {
        return l_p;
    } else {
        return find_most_probable_path(false);
    }
}


bool PolarCode::crc_check(const u8 *info_bit_padded) {
    bool crc_pass = true;
    for (u16 i = info_length; i < info_length + crc_size; ++i) {
        u8 crc_bit = 0;
        for (u16 j = 0; j < info_length; ++j) {
            crc_bit = (u8) ((crc_bit + crc_matrix.at(i - info_length).at(j) *
                                       info_bit_padded[channel_order_descending.at(j)]) % 2);
        }

        if (crc_bit != info_bit_padded[channel_order_descending.at(i)]) {
            crc_pass = false;
            break;
        }
    }

    return crc_pass;
}
