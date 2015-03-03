#include "histo-header-boost.h"
#include "gmock/gmock.h"
#include <memory>
#include <iostream>
#include <random>
using namespace testing;
using namespace std;
using namespace histo;

TEST(HistoBoost, ConstructorBreaks){
    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    vector<double> br{1.0, 2.0, 15.0, 20.0};
    HistoB<double>(data, br);

}
