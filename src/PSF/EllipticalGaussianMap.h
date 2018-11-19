/**\file
 *
 * \brief Declares a smoothly varying over an image elliptical Gaussian PSF
 * map class.
 */

#ifndef __SDK_PSF_MAP_H
#define __SDK_PSF_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EllipticalGaussian.h"
#include "Map.h"
#include "../IO/H5IODataTree.h"
#include "../Core/NaN.h"
#include "Eigen/Dense"

namespace PSF {

    ///Smoothly varying over the image elliptical Gaussian PSF.
    class LIB_PUBLIC EllipticalGaussianMap : public Map {
    private:
        ///Simplify statements requiring this type.
        typedef Eigen::Map<
            const Eigen::Matrix<double,
                                Eigen::Dynamic,
                                Eigen::Dynamic, 
                                Eigen::RowMajor> 
        > RowMajorMap;

        Eigen::MatrixXd
            ///\brief Apply to a vector of polynomial terms to get the values 
            ///of the S, D and K parameters, in that order.
            __sdk_expansion;

    public:
        ///Construct a map from the given data tree.
        EllipticalGaussianMap(
            ///The I/O data tree to construct the map from. It should contain
            ///at least the quantities identified by
            ///EllipticalGaussianMap::required_data_tree_quantities().
            const IO::H5IODataTree &data
        );

        ///Evaluate the map at particular values of the terms.
        EllipticalGaussian *get_psf(
            ///The values of the terms on which the PSF map depends.
            const Eigen::VectorXd &terms,

            ///Background to add to the PSF (assumed constant) normalized in
            ///the same way as the backgroundless PSFs produced by the map.
            double background=0
        ) const;

        ///A reference to a dynamically allocated PSF.
        PSF *operator()(
            ///The values of the terms on which the PSF map depends.
            const Eigen::VectorXd &terms,

            ///Background to add to the PSF (assumed constant) normalized in
            ///the same way as the backgroundless PSFs produced by the map.
            double background=0
        ) const
        {return get_psf(terms, background);}

        ///\brief All quantities needed to construct the PSF map from an I/O 
        ///data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    }; //End EllipticalGaussianMap class.

} //End PSF namespace.

#endif
