/**\file
 *
 * \brief Some useful typedef statements.
 */

#ifndef __TYPEDEFS_H
#define __TYPEDEFS_H

#include <vector>
#include <list>
#include <string>
#include <valarray>
#include "PhotColumns.h"
#include "Eigen/Dense"

namespace Core {

    ///Flags for the quality of a photometric measurement.
    enum PhotometryFlag {UNDEFINED=-1, GOOD, SATURATED, BAD};

    typedef std::vector<double>::size_type vector_size_type;

    ///Synonim for list of doubles (needed for boost command line parsing).
    class RealList : public std::list<double> {};

    ///\brief Synonim for list of column names (needed for boost command 
    ///line parsing).
    class ColumnList : public std::list<Phot::Columns> {};

    ///Synonym for list of strings (needed for boost command line parsing).
    class StringList : public std::list<std::string> {};

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

} //End Core namespace.


#endif
