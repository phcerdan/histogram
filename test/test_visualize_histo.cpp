#include "gmock/gmock.h"
#include "visualize_histo.hpp"
#include <memory>
#include <iostream>
#include <random>
using namespace testing;
using namespace std;
using namespace histo;

TEST(VisualizeHisto, withJustData){

    vector<double> data{0.0, 1.0, 1.0,1.0, 2.0, 3.0, 5.0, 5.0, 8.0, 8.0,  12.0};
    Histo<double> h(data, histo::GenerateBreaksFromRangeAndBins<double>(0.0,15.0, 5));
    h.PrintBreaksAndCounts(std::cout);
    h.name = "withJustData";
    visualize_histo(h);
    visualize_histo(h, vtkChart::BAR);
}
