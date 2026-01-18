#include <gfeat>

#include "gtest/gtest.h"

TEST(TestAOD1B, Read) {
    AOD1B aod1b;
    std::string filename =
        std::string(DATA_DIR) + "/aod1b/AOD1B_2025-01-01_X_07.asc";
    aod1b.load(filename);
    DateTime date(2024, 1, 1, 0);
    auto data = aod1b.get(date, GLO);
    std::cout << data.coefficients << std::endl;
}