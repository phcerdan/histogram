/**
@file histo-header-boost.h
Histogram in C++11, 1D, simple header-only, inspired by R, calculate and optimize this->breaks automatically. Accepts different precissions. Implements Save method using boost::filesystem.
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

#ifndef HISTO_HEADER_BOOST_H_
#define HISTO_HEADER_BOOST_H_
#include "histo-header.h"
#include <fstream>
#include <boost/filesystem.hpp>
/** histo namespace in histo-header.h*/
namespace histo {

template <typename T, typename PRECI = double, typename PRECI_INTEGER = unsigned long int>
struct HistoB : public histo::Histo<T, PRECI, PRECI_INTEGER> {

/********** CONSTRUCTORS ************/
    /// C++11 Constructors inheritance
    using Histo<T, PRECI, PRECI_INTEGER>::Histo;

    void Save(const std::string & fname, const std::string & relativeOutputFolder = "./"){
        boost::filesystem::path dir(relativeOutputFolder);
        boost::filesystem::path histopath(dir/(fname +".histo"));
        if (relativeOutputFolder != "./"){
            boost::filesystem::create_directories(dir);
        }
        std::ofstream output_file (histopath.string()); // delete everything inside the file(default)
        if (!output_file.is_open()) {
            perror( ("Error creating histo file in " + histopath.string() ).c_str());
        }

        for (unsigned long long i = 0; i!=this->counts.size() - 1; i++ ){
            PRECI m_break = (this->breaks[i + 1] - this->breaks[i]) / 2.0;
            output_file << m_break << " " << this->counts[i] << std::endl;
        }

        if(output_file.bad()){
            perror( ("Error saving histo file in " + histopath.string() ).c_str());
        }
        output_file.close();
    }
};

}
#endif
