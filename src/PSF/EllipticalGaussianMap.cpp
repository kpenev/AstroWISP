/**\file
 *
 * \brief Definitions of some methods of the EllipticalGaussianMap class.
 */

#include "EllipticalGaussianMap.h"

namespace PSF {

    EllipticalGaussianMap::EllipticalGaussianMap(
        const IO::H5IODataTree &data
    )
    {
        const std::vector<double>& coefficients =
            data.get< std::vector<double> >(
                    "psffit.psfmap",
                    std::vector<double>(),
                    IO::TranslateToAny< std::vector<double> >()
            );
        assert(coefficients.size() % 3 == 0);
        unsigned nterms = coefficients.size() / 3;
        set_num_terms(nterms);
        __sdk_expansion = RowMajorMap(&(coefficients[0]), 3, nterms);
    }

    EllipticalGaussian *EllipticalGaussianMap::get_psf(
            const Eigen::VectorXd &terms,
            double background
    ) const
    {
        Eigen::VectorXd sdk = __sdk_expansion * terms;
        return new EllipticalGaussian(sdk[0], sdk[1], sdk[2], background);
    }

    const std::set<std::string> &
    EllipticalGaussianMap::required_data_tree_quantities()
    {

        const std::string additional_required_data_tree_quantities[] = {
#ifdef DEBUG
            "psffit.model",
#endif
            "psffit.psfmap"
        };

        static const std::set<std::string> required_data_tree_quantities =
            combine_required_tree_quantities(
                    additional_required_data_tree_quantities,
                    additional_required_data_tree_quantities
                    +
                    sizeof(additional_required_data_tree_quantities)
                    /
                    sizeof(additional_required_data_tree_quantities[0])
            );
        return required_data_tree_quantities;
    }

} //End PSF namespace.
