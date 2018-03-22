#include "gmock/gmock.h"
#include "histo.hpp"
#include <memory>
#include <iostream>
#include <random>
using namespace testing;
using namespace std;
using namespace histo;

TEST(HistoConstructor, withJustData){

    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    Histo<double> h(data); // Default method to calculate Breaks.
    pair<double, double> expected_range = make_pair(1.0, 19.0);
    EXPECT_EQ(expected_range, h.range);
    double sum_counts{0};
    for(auto &c : h.counts){
        sum_counts +=c;
    }
    EXPECT_EQ(5, sum_counts);
    h.PrintBreaksAndCounts(std::cout);
    h.PrintCentersAndCounts(std::cout);
    std::cout << "breaks: ";
    h.PrintBreaks(std::cout);
    std::cout << "counts: ";
    h.PrintCounts(std::cout);
    std::cout << "centers: ";
    h.PrintCenters(std::cout);
}

TEST(HistoConstructor, withInputRange){
    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    std::pair<double,double> input_range = make_pair(-5,24.0);
    Histo<double> h(data, input_range); // Default method to calculate Breaks.
    double sum_counts{0};
    for(auto &c : h.counts){
        sum_counts +=c;
    }
    EXPECT_EQ(5, sum_counts);
}

TEST(HistoConstructor, withBreaks){
    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    vector<double> br{1.0, 2.0, 15.0, 20.0};
    Histo<double> h(data, br);
    EXPECT_EQ(3, h.bins );
    EXPECT_EQ(2, h.counts[0]);
    EXPECT_EQ(2, h.counts[1]);
    EXPECT_EQ(1, h.counts[2]);
    data = {20.0, 20.0, 20.0, 20.0, 20.0};
    h.FillCounts(data);
    EXPECT_EQ(6, h.counts[2]);
    data = {-1.0};
    EXPECT_THROW(h.FillCounts(data),histo_error);
}

TEST(HistoConstructor, withBreaksWithGenerator){
    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    unsigned int input_bins = 10;
    Histo<double> h(data, histo::GenerateBreaksFromRangeAndBins<double>(0.0,20.0, input_bins));
    vector<double> expected_breaks;
    double expected_width = 2.0;
    for(unsigned int i = 0; i!=input_bins + 1; i++){
        expected_breaks.push_back(expected_width * i);
    }
    EXPECT_EQ(expected_breaks, h.breaks);
}

/**
 * @brief Test the constructor with int data and input_breaks.
 */
struct HistoIntInputBreaks : public ::testing::Test{
    static vector<int> data;
    static Histo<double>   h;
};
vector<int> HistoIntInputBreaks::data{-2,-1,0,1,2};
Histo<double> HistoIntInputBreaks::h(data, histo::GenerateBreaksFromRangeAndBins<double>(-2,2,2));

TEST_F(HistoIntInputBreaks, MembersAreCorrect) {
    EXPECT_EQ(2, (int)h.bins);
    std::vector<double> gen_breaks{-2,0,2};
    EXPECT_TRUE(gen_breaks == h.breaks) ;
    EXPECT_TRUE(std::make_pair(-2.0,2.0) == h.range);
    EXPECT_EQ(2, (int)h.counts[0]);
    EXPECT_EQ(3, (int)h.counts[1]);
}
TEST_F(HistoIntInputBreaks, AddMoreDataToExistingHisto) {
    vector<int> extra_data{-1,-1,2};
    h.FillCounts(extra_data);
    EXPECT_EQ(4, (int)h.counts[0]);
    EXPECT_EQ(4, (int)h.counts[1]);
}
/**
 * @brief You can modify counts manually, or with methods Increase, Decrease, SetCount.
 * The latter provides security checks, throwing exepction if exceed limits (big, or less than zero).
 */
TEST_F(HistoIntInputBreaks, ModifyCountsManuallyOrMethods){
    EXPECT_EQ(0, h.IndexFromValue(-2));
    auto old_count = h.counts[0]++;
    EXPECT_EQ(old_count + 1, h.counts[0] );
    // Manipulating counts directly do not throw expections.
    EXPECT_NO_THROW(h.counts[0] = -2);
    EXPECT_NO_THROW(h.counts[0]++);
    // Decrease throws error.
    h.counts[0] = 0;
    EXPECT_THROW(h.Decrease(0), histo_error);
    // Increase throws errors.
    h.counts[0] = numeric_limits<unsigned long int>::max();
    EXPECT_THROW(h.Increase(0), histo_error) << h.counts[0];
    // SetCount throws errors.
    EXPECT_THROW(h.SetCount(100, 2), histo_error) ;
}

/**
 * @brief Test the constructor with just double data and Scott method to generate breaks.
 */
struct HistoDoubleScott : public ::testing::Test{
    static vector<double> data;
    static Histo<double>   h;
};
vector<double> HistoDoubleScott::data{1.0,3.0,5.0};
Histo<double> HistoDoubleScott::h(data, breaks_method::Scott);

TEST_F(HistoDoubleScott, MembersAreCorrect) {
    EXPECT_FLOAT_EQ(1.0, h.range.first);
    EXPECT_FLOAT_EQ(5.0, h.range.second);
    EXPECT_EQ(1, h.bins);
    EXPECT_EQ(1, (int)h.counts.size());
    vector<double> expected_breaks{1.0, 5.0};
    EXPECT_TRUE(expected_breaks == h.breaks);
    EXPECT_EQ(3, h.counts[0]);
    // EXPECT_EQ(2, h.counts[1]);

}
TEST_F(HistoDoubleScott, IndexFromValueWorks) {
    ASSERT_THROW( h.IndexFromValue(1.0 - numeric_limits<double>::epsilon()), histo_error);
    EXPECT_EQ(0, h.IndexFromValue(1.0));
    EXPECT_EQ(0, h.IndexFromValue(h.breaks[1] - 2*numeric_limits<double>::epsilon())) << h.breaks[1];
    EXPECT_EQ(0, h.IndexFromValue(h.breaks[1])) << h.breaks[1];
    ASSERT_THROW(h.IndexFromValue(10), histo_error);
}
/**
 * @brief Test the constructor with long double, input_range, and specific method to generate breaks.
 */
struct HistoLongDoubleInputRange : public ::testing::Test{
    static vector<long double> data;
    static Histo< long double >   h;
};
vector< long double> HistoLongDoubleInputRange::data{-0.8,0.1,0.9};
Histo< long double > HistoLongDoubleInputRange::h(data,std::make_pair(-1.0, 1.0), breaks_method::Scott);

TEST_F(HistoLongDoubleInputRange, BalanceBreaksWorks) {
    std::vector<long double> balanced_breaks{-1.0,1.0};
    EXPECT_EQ(balanced_breaks.size(), h.breaks.size());
    EXPECT_TRUE(balanced_breaks == h.breaks);
}

/**
 * @brief LotsOfData
 */
static default_random_engine generator;
static uniform_real_distribution<double> cosined(-1.0,1.0);

TEST(HistoLotsOfDataInputRange, BalanceBreaksWorks) {
    unsigned long int ndata = 10000;
    std::vector<long double> data(ndata);
    for(auto &x : data){
        x = cosined(generator);
    }
    std::pair<long double, long double> input_pair = std::make_pair(-1.0, 1.0);
    Histo< long double, unsigned long long> h(data, input_pair, breaks_method::Scott);
    EXPECT_FLOAT_EQ(input_pair.first , h.breaks[0]);
    EXPECT_FLOAT_EQ(input_pair.second , h.breaks[h.bins]);
}

TEST(GenerateBreaksFromRangeAndWidth, withSameUpper) {
    double low = 0.0;
    double upper = 4.0;
    double width = 1.0;
    auto breaks = histo::GenerateBreaksFromRangeAndWidth<double>(low, upper, width);
    EXPECT_EQ(breaks.size(), 5);
    EXPECT_FLOAT_EQ(breaks[0], low);
    EXPECT_FLOAT_EQ(breaks[1], 1.0);
    EXPECT_FLOAT_EQ(breaks[2], 2.0);
    EXPECT_FLOAT_EQ(breaks[3], 3.0);
    EXPECT_FLOAT_EQ(breaks[4], 4.0);
}

TEST(GenerateBreaksFromRangeAndWidth, withGreaterUpper) {
    double low = 0.0;
    double upper = 4.5;
    double width = 1.0;
    auto breaks = histo::GenerateBreaksFromRangeAndWidth<double>(low, upper, width);
    EXPECT_EQ(breaks.size(), 6);
    EXPECT_FLOAT_EQ(breaks[0], low);
    EXPECT_FLOAT_EQ(breaks[1], 1.0);
    EXPECT_FLOAT_EQ(breaks[2], 2.0);
    EXPECT_FLOAT_EQ(breaks[3], 3.0);
    EXPECT_FLOAT_EQ(breaks[4], 4.0);
    EXPECT_FLOAT_EQ(breaks[5], 5.0);
}

TEST(NormalizeByArea, withJustData ) {
    double low = 0.0;
    double upper = 20.0;
    double width = 1.0;
    auto breaks = histo::GenerateBreaksFromRangeAndWidth<double>(low, upper, width);
    EXPECT_EQ(breaks.size(), 21);
    vector<double> data{1.0,1.0,2.0, 3.0, 19.0};
    Histo<double> h(data, breaks);
    auto h_norm = NormalizeByArea(h);
    const auto & counts = h_norm.counts;
    double sum_areas = 2*width + (1 * width) * 3;
    EXPECT_FLOAT_EQ(counts[0], 0.0);
    EXPECT_FLOAT_EQ(counts[1], 2.0 / sum_areas);
    EXPECT_FLOAT_EQ(counts[2], 1.0 / sum_areas);
    EXPECT_FLOAT_EQ(counts[3], 1.0 / sum_areas);
    EXPECT_FLOAT_EQ(counts[19], 1.0/ sum_areas);
}
