/**
@file histo-header.h
Histogram in C++11, 1D, simple header-only, inspired by R, calculate and optimize breaks automatically. Accepts different precissions.
Copyright (C) 2015 Pablo Hernandez. github/phcerdan/histo-header.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef HISTO_HEADER_H_
#define HISTO_HEADER_H_
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <iostream>
#ifdef WITH_BOOST
#include <fstream>
#include <boost/filesystem.hpp>
#endif
/** histo namespace in histo-header.h*/
namespace histo {
/** \defgroup breaks_methods breaks_methods */
/**@{
 * @brief Breaks method to optimal calculation of breaks based on input data and range.
 * Scott = 0.
 */
enum breaks_method {
    /** Scott Method*/
    Scott = 0
};
/** @} */

/** \defgroup GenerateBreaks Generate breaks from data, range, and/or bins. */
/** @{
 * @brief Help functions to manually creating breaks from input range
 * (low, upper) and desired number of bins.
 *
 * @tparam T Input data type.
 * @tparam PRECI It should be equal to T, except when T is int.
 * @param low first value of breaks.
 * @param upper last value of breaks.
 * @param bins number of divisions.
 *
 * @return breaks vector with the frontier values of each bin.
 */
template<typename T, typename PRECI = double>
std::vector<T> GenerateBreaksFromRangeAndBins(const T& low, const T& upper, const unsigned long int &bins){
    std::vector<T> breaks(bins+1);
    T width = (upper - low)/ static_cast<PRECI>(bins);
    for (unsigned long int i = 0; i!=bins+1; i++){
        breaks[i] = low + i*width;
    }
    return breaks;
};

/** @brief @sa GenerateBreaksFromRangeAndBins() */
template<typename T, typename PRECI = double>
std::vector<T> GenerateBreaksFromRangeAndBins(const std::pair<T,T> &range_low_upper, const unsigned long int &bins){
   T low   = range_low_upper.first;
   T upper = range_low_upper.second;
   return GenerateBreaksFromRangeAndBins<T, PRECI>(low, upper, bins);
};
/** @} */

/** Exception class for  Histo */
class histo_error : public std::runtime_error {
public:
    /** Error constructor from runtime_error*/
    histo_error(const std::string & s) : std::runtime_error(s){ };
};


/** @brief Variance calculation from Container with data.
 *
 * @tparam T Type of data.
 * @tparam Container std:: type containing data (vector, array,...)
 * @param xs The container.
 *
 * @return Variance of type T.
 */
template <typename T, typename Container>
T variance_welford(const Container& xs)
{
    unsigned long long N = 0;
    T M = 0, S = 0, Mprev = 0;
    for(auto x : xs) {
        ++N;
        Mprev = M;
        M += (x - Mprev) / N;
        S += (x - Mprev) * (x - M);
    }
    return S / (N-1);
}

/** Precise comparison: equal than.
* @param v1 type T, variable 1 to compare
* @param v2 type T, variable 2 to compare
* @return bool
*/
template<typename T, unsigned int N = 1>
bool isequalthan(const T& v1, const T& v2)
{
    return std::abs(v1-v2)<= N * std::numeric_limits<T>::epsilon();
}


/**
 * @brief Histogram inspired by R.
 * Simple, no dependancies, header-only.
 * It accepts different data values and precission.
 *
 * @tparam T Input data type.
 * @tparam PRECI It should be equal to T, except when T is int.
 * @tparam PRECI_INTEGER int type for counts data member.
 */
template <typename T, typename PRECI = double, typename PRECI_INTEGER = unsigned long int>
struct Histo {

/************* DATA *****************/
    /** Low and upper limit for breaks. */
    std::pair<T,T> range;
    /** Value of the breaks between bins, [low,...,upper].
     *  size.breaks = size.counts + 1. */
    std::vector<T> breaks;
    /** breaks.size() */
    unsigned long int bins{0};
    /** int vector holding the counts for each breaks interval.*/
    std::vector<PRECI_INTEGER> counts;

/********** CONSTRUCTORS ************/
    Histo() = default;

    /**
     * @brief Constructor that takes range as the min, and max values of data.
     *
     * @param data
     * @param method Method to calculate breaks from @sa histo::breaks_method
     */
    Histo( const std::vector<T> &data, histo::breaks_method method = Scott )
    {
        auto range_ptr = std::minmax_element(data.begin(), data.end());
        range     = std::make_pair(*range_ptr.first, *range_ptr.second);
        breaks    = CalculateBreaks(data, range, method);
        bins      = static_cast<decltype(bins)>( breaks.size() - 1);
        ResetCounts();
        FillCounts(data);
    };

    /**
     * @brief Constructor with fixed input range.
     * @param data
     * @param input_range low and upper value
     * @param method Method to calculate breaks from @sa histo::breaks_method
     */
    Histo(const std::vector<T> &data, const std::pair<T,T> &input_range, histo::breaks_method method = Scott ){
        range     = input_range;
        breaks    = CalculateBreaks(data, range, method);
        bins      = static_cast<decltype(bins)>( breaks.size() - 1 );
        ResetCounts();
        FillCounts(data);
    };
    /**
     * @brief Constructor that accepts a vector of breaks.
     * You can use @sa histo::GenerateBreaksFromRangeAndBins
     * to help you creating the vector from specific number of bins and range.
     *
     * @param data
     * @param input_breaks
     */
    Histo(const std::vector<T> &data, const std::vector<T> &input_breaks) {
        breaks = input_breaks;
        if(!CheckIfMonotonicallyIncreasing(breaks)) throw
            histo_error("input_breaks are not monotocally increasing");
        range  = std::make_pair(breaks[0], breaks[breaks.size() - 1]);
        bins    = static_cast<decltype(bins)>( breaks.size() - 1 );
        ResetCounts();
        FillCounts(data);
    };

/********* PUBLIC METHODS ***********/
    /**
     * @brief Return the index of @sa counts associated to the input value
     *
     * @param value Ranging from range.first to range.second
     * @return Index of counts
     */
    unsigned long int IndexFromValue(const T &value){
        std::vector<float>::iterator low = std::lower_bound(breaks.begin(), breaks.end(), value);
        return low - breaks.begin();
        //unsigned long int lo{0},hi{bins}, newb; // include right border in the last bin.
        //if(value >= breaks[lo] && (value < breaks[hi] || histo::isequalthan<T>(value,breaks[hi]) )){
        //    while( hi - lo >= 2){
        //        newb = (hi+lo)/2;
        //        if ( (value >= breaks[newb]) ) lo = newb;
        //        else hi = newb;
        //    }
        //} else {
        //    throw histo_error(" IndexFromValue: "+ std::to_string(value) +  " is out of bonds");
        //}
        //
        //return lo;
    };

    /** @brief Resize counts and reset value to zero. */
    void ResetCounts(){
        counts.resize(bins);
        for (auto &c : counts){
            c = 0;
        }
    };
    /**
     * @brief Fill counts from data.
     * Breaks must have been set-up before calling this method.
     *
     * @param data
     *
     * @return Reference to the data member @sa counts
     */
    std::vector<PRECI_INTEGER>& FillCounts(const std::vector<T> &data){
        for (auto &v : data){
            counts[IndexFromValue(v)] ++ ;
        }
        return counts;
    };

    /** \defgroup CountsManipulation Counts Safe Manipulation */
    /** @{
     * @brief Increase count by one, checking if exceeds max_integer_.
     * @param index of counts
     */
    void Increase(const unsigned long int & index){
        if (counts[index] == max_integer_)
            throw histo_error("Increase has exceded PRECI_INTEGER."
                    " Index: " + std::to_string(index) + " Value: " + std::to_string(counts[index]) );
        counts[index]++;
    };

    /** @brief Decrease count by one, checking if it goes negative.
     * @param index of counts.
     */
    void Decrease(const unsigned long int & index){
        if (counts[index] <= 0 )
            throw histo_error("Decrease has reached negative value."
                    " Index: " + std::to_string(index) + " Value: " + std::to_string(counts[index]) );
        counts[index]--;
    };

    /** @brief Set count value. Checks for negative or greater than max_integer_.
     * @param index of counts.
     * @param v value to set.
     */
    void SetCount(const unsigned long int & index, const long double & v ){
        if (v < 0 || v > max_integer_)
            throw histo_error("SetCount to a negative value, or greater than allowed"
                    " Index: " + std::to_string(index) + " Value: " + std::to_string(counts[index]) );
        counts[index] = static_cast<PRECI_INTEGER>(v);
    };

    /** @} */
protected:
    /** Max integer for the current PRECI_INTEGER type*/
    PRECI_INTEGER max_integer_{std::numeric_limits<PRECI_INTEGER>::max()};

    bool CheckIfMonotonicallyIncreasing(const std::vector<T> &input_breaks){
        T prev_value = input_breaks[0];
        for( auto it = input_breaks.begin() + 1, it_end = input_breaks.end();
                it!=it_end; it++)
        {
            if(*it <= prev_value) return false;
            prev_value = *it;
        }
        return true;
    }
    /**
     * @brief Method to wrap breaks calculation methods that
     * take into account the input data and range to optimize breaks vector.
     *
     * @param data
     * @param rang Range of breaks vector (low, upper)
     * @param method Method to calculate breaks from @histo::breaks_method
     *
     * @return Reference to data member: breaks.
     */
    std::vector<T>& CalculateBreaks(const std::vector<T> & data, const std::pair<T,T> & rang, histo::breaks_method method ){
        switch(method) {
        case Scott:
             return ScottMethod(data,rang);
             break;
        default:
             throw histo_error("CalculateBreaks: No Valid Method selected to calculate breaks.");
        }
    };

    bool CheckBreaksAreEquidistant(const std::vector<T> & input_breaks){
        T diff = input_breaks[1] - input_breaks[0];
        for( auto it = input_breaks.begin() + 1, it_end = input_breaks.end();
                it!=it_end; it++)
        {
            // T new_diff = *it - *(it -1);
            // std::cout<<  new_diff <<" " << diff << std::endl;
            // Soft comparisson, high number of epsilons.
            if(!isequalthan<T, 25>(*it - *(it-1), diff)) return false;
        }
        return true;
    };

    bool BalanceBreaksWithRange(std::vector<T> &input_breaks, std::pair<T,T> input_range){
        if (!CheckBreaksAreEquidistant(input_breaks))
            throw histo_error("BalanceBreaksWithRange cannot be applied in NON Equidistant breaks");

        unsigned long int nbins = input_breaks.size() - 1;
        T width                 = input_breaks[1] - input_breaks[0];
        // diff_low is > 0 when it does not reach range, and < 0 when it goes beyond.
        T diff_low              = input_breaks[0] - input_range.first;
        // diff_upper is < 0 when it does not reach range, and >0 when it goes beyond.
        T diff_upper            = input_breaks[nbins] - input_range.second;
        bool diff_low_isZero   = isequalthan<T>(diff_low, 0);
        bool diff_upper_isZero = isequalthan<T>(diff_upper, 0);

        if (diff_low_isZero && diff_upper_isZero) return false;

        // Put the first break in the range.first, and move all breaks accordingly.
        if(!diff_low_isZero) ShiftBreaks(input_breaks, -diff_low);
        CheckAndUpdateDiff(diff_upper, diff_upper_isZero, input_breaks[nbins], input_range.second);
        if(diff_upper_isZero)  return true;

        // Check if it is better to remove the last break than to add more.
        // Biased to add more (1 is no bias).
        double bias_to_add_bin = 0.8;
        T diff_upper_before = input_breaks[nbins - 1] - input_range.second;
        if ((diff_upper_before < 0) && (diff_upper > 0 )
                && ( fabs(diff_upper_before) < bias_to_add_bin * fabs(diff_upper) ))
        {
            nbins--;
            input_breaks.pop_back();
            T width_to_expand = diff_upper_before/nbins;
            ShrinkOrExpandBreaks(input_breaks, -width_to_expand);
        }
        CheckAndUpdateDiff(diff_upper, diff_upper_isZero, input_breaks[nbins], input_range.second);
        if(diff_upper_isZero)  return true;

        // If diff_upper is <0,  add bins until reach the upper range.
        while(diff_upper < 0){
            nbins++;
            input_breaks.push_back(input_range.first + nbins * width);
            diff_upper  = input_breaks[nbins] - input_range.second;
        }
        diff_upper_isZero = isequalthan<T>(diff_upper, 0);
        if (diff_upper_isZero) return true;

        // If diff_upper > 0.
        // Shrink the width of the breaks, and calculate the new positions.
        T width_to_shrink = diff_upper / nbins;
        ShrinkOrExpandBreaks(input_breaks, - width_to_shrink);
        return true;
    };

    void CheckAndUpdateDiff(T & diff, bool & diff_isZero, const T & rhs, const T & lhs){
        diff = rhs - lhs;
        diff_isZero = isequalthan<T>(diff, 0);
    };
    void ShiftBreaks(std::vector<T> &input_breaks, const T & d){
        for (auto &v : input_breaks){
            v = v + d;
        }
    };
    void ShrinkOrExpandBreaks(std::vector<T> &input_breaks, const T & d){
        unsigned long int i{0} ;
        for (auto &v : input_breaks){
            v = v + i*d;
            i++;
        }
    };
    /**
     * @brief Scott Method to calculate optimal breaks.
     * Calculate variance from @histo::variance_welford
     *
     * @param data
     * @param rang Range of breaks vector (low, upper)
     *
     * @return Reference to data member: breaks.
     */
    std::vector<T>& ScottMethod(const std::vector<T> &data, const std::pair<T,T> &rang){
        PRECI sigma = variance_welford<PRECI,std::vector<T>>(data);
        PRECI width  = 3.5 * sqrt(sigma) / static_cast<PRECI>( data.size() );
        bins    = std::ceil( (rang.second - rang.first) / width);
        breaks.resize(bins + 1 );
        for (unsigned long int i = 0; i!=bins + 1; i++){
            breaks[i] = range.first + i*width;
        }
        // std::cout << "Non balanced breaks" << std::endl;
        // std::cout<< "bins is: " << bins <<" width is: "<< width << " sigma is: " << sigma  <<std::endl;
        // std::cout << "first break" << breaks[0] << " Second Break" << breaks[bins]  << std::endl;
        // std::for_each(std::begin(breaks), std::end(breaks), [](const T& v) {std::cout<<v <<std::endl;});

        BalanceBreaksWithRange(breaks, range);
        bins = breaks.size() - 1;
        // std::cout << "Balanced new breaks" << std::endl;
        // std::cout<< "bins is: " << bins <<" width is: "<< breaks[1]-breaks[0]  <<std::endl;
        // std::cout << "first break" << breaks[0] << " Second Break" << breaks[bins]  << std::endl;
        // std::for_each(std::begin(breaks), std::end(breaks), [](const T& v) {std::cout<<v <<std::endl;});
        return breaks;
    };
};

} // End of namespace histo
#endif