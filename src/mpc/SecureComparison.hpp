#ifndef MPC_SECURE_COMPARISON_HPP
#define MPC_SECURE_COMPARISON_HPP

#include <vector>
#include <iomanip>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <limits>

void SecureComparison(int a_silo_id, int b_silo_id, float a_value, float b_value, int& silo_id);

void SecureComparison(int a_silo_id, int b_silo_id, float a_value, float b_value, int& silo_id) {
    silo_id = (a_value <= b_value) ? a_silo_id : b_silo_id;
}

#endif // MPC_SECURE_COMPARISON_HPP