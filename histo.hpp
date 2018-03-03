/**
@file histo.hpp
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

#ifndef HISTO_HPP_
#define HISTO_HPP_
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <iostream>
/** histo namespace in histo.h*/
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
 * You can use this function to force breaks in the place you want.
 * It is specially useful if you want to store integers.
 * For example, image your data is unsigned, [0,...,10]
 * You can choose breaks in a way that the middle of the bins are unsigned:
 * GenerateBreaksFromRangeAndBins( min - 0.5, max + 0.5, max + 1)
 * If you use PrintCenter, you will the expected histo.
 *
 * @tparam PRECI is the type of the breaks, should be greater or as precise as T.
 * @param low first value of breaks.
 * @param upper last value of breaks.
 * @param bins number of divisions.
 *
 * @return breaks vector with the frontier values of each bin.
 */
template<typename PRECI = double>
std::vector<PRECI> GenerateBreaksFromRangeAndBins(
      const PRECI& low, const PRECI& upper, const unsigned long int &bins){
    std::vector<PRECI> breaks(bins+1);
    PRECI width = (upper - low)/ static_cast<PRECI>(bins);
    for (unsigned long int i = 0; i!=bins+1; i++){
        breaks[i] = low + i*width;
    }
    return breaks;
};

/** @brief @sa GenerateBreaksFromRangeAndBins() */
template<typename PRECI = double>
std::vector<PRECI> GenerateBreaksFromRangeAndBins(
      const std::pair<PRECI,PRECI> &range_low_upper, const unsigned long int &bins){
   auto low   = range_low_upper.first;
   auto upper = range_low_upper.second;
   return GenerateBreaksFromRangeAndBins<PRECI>(low, upper, bins);
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
 * @tparam TData Type of data.
 * @tparam Container std:: type containing data (vector, array,...)
 * @param xs The container.
 *
 * @return Variance of type TData
 */
template <typename TData, typename Container>
TData variance_welford(const Container& xs)
{
    unsigned long long N = 0;
    TData M = 0, S = 0, Mprev = 0;
    for(const auto & x : xs) {
        ++N;
        Mprev = M;
        M += (x - Mprev) / N;
        S += (x - Mprev) * (x - M);
    }
    return S / (N-1);
}

/** Precise comparison: equal than.
* @param v1 type TData, variable 1 to compare
* @param v2 type TData, variable 2 to compare
* @return bool
*/
template<typename TData, unsigned int N = 1>
bool isequalthan(const TData& v1, const TData& v2)
{
    return std::abs(v1-v2)<= N * std::numeric_limits<TData>::epsilon();
}


/**
 * @brief Histogram inspired by R.
 * Simple, no dependancies, header-only.
 * It accepts different data values and precission.
 *
//  * @tparam T Input data type.
 * @tparam PRECI It should be equal to T, except when T is int.
 * @tparam PRECI_INTEGER int type for counts data member.
 */
template <typename PRECI = double, typename PRECI_INTEGER = unsigned long int>
struct Histo {

/************* DATA *****************/
    /** Low and upper limit for breaks. */
    std::pair<PRECI,PRECI> range;
    /** Value of the breaks between bins, [low,...,upper].
     *  size.breaks = size.counts + 1. */
    std::vector<PRECI> breaks;
    /** breaks.size() - 1 */
    unsigned long int bins{0};
    /** int vector holding the counts for each breaks interval.*/
    std::vector<PRECI_INTEGER> counts;
    /** name/description of the histogram */
    std::string name;

/********** CONSTRUCTORS ************/
    Histo() = default;

    /**
     * @brief Constructor that takes range as the min, and max values of data.
     *
     * @param data
     * @param method Method to calculate breaks from @sa histo::breaks_method
     */
    template<typename TData>
    Histo( const std::vector<TData> &data, histo::breaks_method method = Scott )
    {
        auto range_ptr = std::minmax_element(data.begin(), data.end());
        range     = std::make_pair( static_cast<PRECI>(*range_ptr.first),
                                    static_cast<PRECI>(*range_ptr.second));
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
    template<typename TData>
    Histo(const std::vector<TData> &data, const std::pair<PRECI,PRECI> &input_range,
          histo::breaks_method method = Scott ){
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
    template<typename TData>
    Histo(const std::vector<TData> &data, const std::vector<PRECI> &input_breaks) {
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
     * @brief print to input std::ostream breaks and counts
     *
     * @param os input ostream, std::cout, std::ofstream, etc.
     */
    void PrintBreaksAndCounts( std::ostream & os ){
       os.setf(std::ios_base::fixed, std::ios_base::floatfield);
       os.precision(9);
       for (unsigned long long i = 0; i < this->counts.size(); i++ ){
          os << "[";
          os << std::setw(18) << this->breaks[i] << "," <<
                std::setw(18) << this->breaks[i+1];
          if (i == counts.size() - 1)
             os << "]";
          else
             os << ")";
          os << " " << std::setw(18) << this->counts[i] << std::endl;
       }
    }

    /**
     * @brief print to input std::ostream center of bins and counts
     *
     * @param os input ostream, std::cout, std::ofstream, etc.
     */
    void PrintCentersAndCounts( std::ostream & os ){
       os.setf(std::ios_base::fixed, std::ios_base::floatfield);
       os.precision(9);
       for (unsigned long long i = 0; i < this->counts.size(); i++ ){
          double break_width = (this->breaks[i + 1] - this->breaks[i]) / 2.0;
          os << std::setw(18) << this->breaks[i] + break_width <<
             " " << std::setw(18) << this->counts[i] << std::endl;
       }
    }

    void PrintCenters( std::ostream & os ){
       os.setf(std::ios_base::fixed, std::ios_base::floatfield);
       os.precision(9);
       for (size_t i = 0; i < this->counts.size(); i++ ){
          double break_width = (this->breaks[i + 1] - this->breaks[i]) / 2.0;
          os << std::setw(18) << this->breaks[i] + break_width;
          if (i != this->counts.size() - 1)
             os << " ";
       }
       os << std::endl;
    }

    void PrintBreaks( std::ostream & os ){
       for (size_t i = 0; i < this->breaks.size(); i++ ){
          os << std::setw(18) << std::setprecision(9) << this->breaks[i];
          if (i != this->breaks.size() - 1)
             os << " ";
       }
       os << std::endl;
    }

    void PrintCounts( std::ostream & os ){
       for (size_t i = 0; i < this->counts.size(); i++ ){
          os << std::setw(18) << this->counts[i];
          if (i != this->counts.size() - 1)
             os << " ";
       }
       os << std::endl;
    }
    /**
     * @brief Return the index of @sa counts associated to the input value
     *
     * @param value Ranging from range.first to range.second
     * @return Index of counts
     */
    template<typename TData>
    unsigned long int IndexFromValue(const TData &value){
       // We could use this with a custom comparator:
       // typename std::vector<T>::iterator low = std::lower_bound(breaks.begin(), breaks.end(), value);
        unsigned long int lo{0},hi{bins}, newb; // include right border in the last bin.
        if(value >= breaks[lo] && (value < breaks[hi] ||
                 histo::isequalthan<PRECI>(value,breaks[hi]) )){
           while( hi - lo >= 2){
               newb = (hi+lo)/2;
               if ( (value >= breaks[newb]) ) lo = newb;
               else hi = newb;
           }
        } else {
           throw histo_error(" IndexFromValue: "+ std::to_string(value) +  " is out of bonds");
        }

        return lo;
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
    template<typename TData>
    std::vector<PRECI_INTEGER>& FillCounts(const std::vector<TData> &data){
        for (auto &v : data){
            counts[IndexFromValue(v)] ++ ;
        }
        return counts;
    };

    /** \defgroup CountsManipulation Counts Safe Manipulation */
    /** @{
     * @brief Increase count by one, checking if exceeds
     * std::numeric_limits<PRECI_INTEGER>::max().
     *
     * @param index of counts
     */
    void Increase(const unsigned long int & index){
        if (counts[index] == std::numeric_limits<PRECI_INTEGER>::max())
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

    /** @brief Set count value. Checks for index being out of bounds.
     * @param index of counts.
     * @param v value to set.
     */
    void SetCount(const unsigned long int & index, const PRECI_INTEGER & v ){
        if (index > bins)
            throw histo_error("Index is out of bounds in SetCount"
                    " Index: " + std::to_string(index) + " Max Bins: " + std::to_string(bins) );
        counts[index] = v;
    };

    /** @} */
protected:
    bool CheckIfMonotonicallyIncreasing(const std::vector<PRECI> &input_breaks){
        auto prev_value = input_breaks[0];
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
    template<typename TData>
    std::vector<PRECI>& CalculateBreaks(const std::vector<TData> & data,
          const std::pair<PRECI,PRECI> & rang, histo::breaks_method method ){
        switch(method) {
        case Scott:
             return ScottMethod(data,rang);
             break;
        default:
             throw histo_error("CalculateBreaks: No Valid Method selected to calculate breaks.");
        }
    };

    bool CheckBreaksAreEquidistant(const std::vector<PRECI> & input_breaks){
        PRECI diff = input_breaks[1] - input_breaks[0];
        for( auto it = input_breaks.begin() + 1, it_end = input_breaks.end();
                it!=it_end; it++)
        {
            // auto new_diff = *it - *(it -1);
            // std::cout<<  new_diff <<" " << diff << std::endl;
            // Soft comparisson, high number of epsilons.
            if(!isequalthan<PRECI, 25>(*it - *(it-1), diff)) return false;
        }
        return true;
    };

    bool BalanceBreaksWithRange(std::vector<PRECI> &input_breaks, std::pair<PRECI,PRECI> input_range){
        if (!CheckBreaksAreEquidistant(input_breaks)){
           std::ostream_iterator<PRECI> out_it (std::cerr,", ");
           std::copy ( input_breaks.begin(), input_breaks.end(), out_it );
           throw histo_error("BalanceBreaksWithRange cannot be applied in NON Equidistant breaks");
        }

        unsigned long int nbins = input_breaks.size() - 1;
        PRECI width                 = input_breaks[1] - input_breaks[0];
        // diff_low is > 0 when it does not reach range, and < 0 when it goes beyond.
        PRECI diff_low              = input_breaks[0] - input_range.first;
        // diff_upper is < 0 when it does not reach range, and >0 when it goes beyond.
        PRECI diff_upper            = input_breaks[nbins] - input_range.second;
        bool diff_low_isZero   = isequalthan<PRECI>(diff_low, 0);
        bool diff_upper_isZero = isequalthan<PRECI>(diff_upper, 0);

        if (diff_low_isZero && diff_upper_isZero) return false;

        // Put the first break in the range.first, and move all breaks accordingly.
        if(!diff_low_isZero) ShiftBreaks(input_breaks, -diff_low);
        CheckAndUpdateDiff(diff_upper, diff_upper_isZero, input_breaks[nbins], input_range.second);
        if(diff_upper_isZero)  return true;

        // Check if it is better to remove the last break than to add more.
        // Biased to add more (1 is no bias).
        double bias_to_add_bin = 0.8;
        PRECI diff_upper_before = input_breaks[nbins - 1] - input_range.second;
        if ((diff_upper_before < 0) && (diff_upper > 0 )
                && ( fabs(diff_upper_before) < bias_to_add_bin * fabs(diff_upper) ))
        {
            nbins--;
            input_breaks.pop_back();
            PRECI width_to_expand = diff_upper_before/nbins;
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
        diff_upper_isZero = isequalthan<PRECI>(diff_upper, 0);
        if (diff_upper_isZero) return true;

        // If diff_upper > 0.
        // Shrink the width of the breaks, and calculate the new positions.
        PRECI width_to_shrink = diff_upper / nbins;
        ShrinkOrExpandBreaks(input_breaks, - width_to_shrink);
        return true;
    };

    void CheckAndUpdateDiff(PRECI & diff, bool & diff_isZero, const PRECI & rhs,
          const PRECI & lhs){
        diff = rhs - lhs;
        diff_isZero = isequalthan<PRECI>(diff, 0);
    };
    void ShiftBreaks(std::vector<PRECI> &input_breaks, const PRECI & d){
        for (auto &v : input_breaks){
            v = v + d;
        }
    };
    void ShrinkOrExpandBreaks(std::vector<PRECI> &input_breaks, const PRECI & d){
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
    template<typename TData>
    std::vector<PRECI>& ScottMethod(const std::vector<TData> &data,
          const std::pair<PRECI,PRECI> &rang){
        PRECI sigma = variance_welford<PRECI>(data);
        // cbrt is cubic root
        PRECI width  = 3.5 * sqrt(sigma) / std::cbrt(static_cast<PRECI>( data.size() ) );
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
