#include <gfeat>

#include "gtest/gtest.h"

TEST(TestAOD1B, Read) {
    AOD1B aod1b;
    aod1b.load("aod1b/AOD1B_2025-01-01_X_07.asc");
    DateTime date(2025, 1, 1, 0);
    auto data = aod1b.get(date, ATM);
}