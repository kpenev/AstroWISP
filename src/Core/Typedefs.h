/**\file
 *
 * \brief Some useful typedef statements.
 */

#ifndef __TYPEDEFS_H
#define __TYPEDEFS_H

#include "../Core/SharedLibraryExportMacros.h"
#include <vector>
#include <list>
#include <string>
#include <valarray>
#include "PhotColumns.h"
#include "Eigen/Dense"

namespace Core {

    ///Flags for the quality of a photometric measurement.
    enum LIB_PUBLIC PhotometryFlag {UNDEFINED=-1, GOOD, SATURATED, BAD};

    ///Alias for the type to use for size/index of double vectors.
    typedef std::vector<double>::size_type vector_size_type;

    ///Synonim for list of doubles (needed for boost command line parsing).
    class LIB_LOCAL RealList : public std::list<double> {};

    ///\brief Synonim for list of column names (needed for boost command
    ///line parsing).
    class LIB_LOCAL ColumnList : public std::list<Phot::Columns> {};

    ///Synonym for list of strings (needed for boost command line parsing).
    class LIB_LOCAL StringList : public std::list<std::string> {};

    ///An Eigen integer matrix suitable organized for saving as FITS image.
    typedef Eigen::Matrix<int,
                          Eigen::Dynamic,
                          Eigen::Dynamic,
                          Eigen::ColMajor> IntImageMatrix;

    ///An Eigen double matrix suitable organized for saving as FITS image.
    typedef Eigen::Matrix<double,
                          Eigen::Dynamic,
                          Eigen::Dynamic,
                          Eigen::ColMajor> DoubleImageMatrix;
    ///
    class DoubleValarray : public std::valarray<double> {
    public:
      DoubleValarray & operator>>=(const double&){return *this;}
      DoubleValarray & operator>>=(const DoubleValarray&){return *this;}
    };

    ///When valarray is done, it tries to do two versions of the operator, the one that takes the entire array is defined in std but it apparently is falling back on the valarray version of it, so maybe change it to take std::valarray<double>&
} //End Core namespace.


#endif
