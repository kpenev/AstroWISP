/** \file
 *
 * \brief The executable that corrects for sub-pixel sensitivity variation.
 *
 * The SubPixPhot executable performs aperture photometry on a list of
 * sources with known point spread functions. It properly accounts for pixels
 * that are only partially inside the aperture as well as non-uniform
 * sensitivity over a pixel (if a sub-pixel sensitivity map is provided).
 *
 * \ingroup SubPixPhot
 */

#include "SubPixPhotIO.h"
#include "../PSF/EllipticalGaussianMap.h"
#include "../PSF/PiecewiseBicubicMap.h"
#include "../PSF/DataTreeCalculations.h"
#include "../PSF/CommandLineUtil.h"
#include "../Background/MeasureAnnulus.h"
#include "../Background/Annulus.h"
#include "../Background/CommandLineUtil.h"
#include "../IO/CommandLineConfig.h"
#include "../IO/FitsImage.h"
#include "../Core/SubPixelCorrectedFlux.h"
#include "../Core/SubPixelMap.h"
#include "../Core/Source.h"
#include "../Core/PhotColumns.h"
#include "../Core/Error.h"
#include "../Core/Typedefs.h"
#include "../Core/CommandLineUtil.h"
#include "Config.h"
#include "Common.h"
#include <string>
#include <list>
#include <ctime>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>


///The hash identifying a unique version of the file.
const std::string SUB_PIX_PHOT_VERSION="$Id: $";

namespace SubPixPhot {

    ///\brief Measures the background of all input sources and adds it to an
    ///I/O tree.
    void add_background_measurements(
        ///The annulus to use for measuring the background.
        const Background::Annulus &annulus,

        ///The image to process
        IO::FitsImage<double> &image,

        ///The I/O tree which must contain the source positions and which
        ///will be updated with the background values, errors etc.
        IO::H5IODataTree &data_tree
    )
    {
        Background::MeasureAnnulus measure_background(
            annulus.inner_radius(),
            annulus.outer_radius(),
            annulus.inner_radius(),
            image
        );
        IO::OutputArray<double> x(data_tree.get<boost::any>("projsrc.x")),
                                y(data_tree.get<boost::any>("projsrc.y"));
        unsigned num_sources = x.size();
        assert(y.size() == num_sources);
        std::vector<double>
            *background = new std::vector<double>(num_sources),
            *background_error = new std::vector<double>(num_sources);
        std::vector<unsigned>
            *background_pixels = new std::vector<unsigned>(num_sources);
        for(
            unsigned source_index = 0;
            source_index < num_sources;
            ++source_index
        )
            measure_background.add_source(x[source_index], y[source_index]);
        measure_background.jump_to_first_source();
        for(
            unsigned source_index = 0;
            source_index < num_sources;
            ++source_index
        ) {
            const Background::Source &this_background = measure_background();
            (*background)[source_index] = this_background.value();
            (*background_error)[source_index] = this_background.error();
            (*background_pixels)[source_index] = this_background.pixels();
            measure_background.next_source();
        }

        IO::TranslateToAny< std::vector<double> > double_trans;
        IO::TranslateToAny< std::vector<unsigned> > unsigned_trans;
        data_tree.put("bg.value", background, double_trans);
        data_tree.put("bg.error", background_error, double_trans);
        data_tree.put("bg.npix", background_pixels, unsigned_trans);
    }

    ///Returns a newly allocated PSF map to use.
    PSF::Map *get_psf_map(
            ///The name of the file to read the PSF map from.
            const std::string &psfmap_filename,

            ///The largest aperture that will be used.
            double max_aperture,

            ///The data tree to fill with extra quantities
            IO::H5IODataTree &psf_data)
    {
        std::set<std::string> quantities_to_read(
            PSF::PiecewiseBicubicMap::required_data_tree_quantities()
        );
        quantities_to_read.insert(
            PSF::EllipticalGaussianMap::required_data_tree_quantities().begin(),
            PSF::EllipticalGaussianMap::required_data_tree_quantities().end()
        );
        quantities_to_read.insert(
            PSF::MapSourceContainer::required_data_tree_quantities().begin(),
            PSF::MapSourceContainer::required_data_tree_quantities().end()
        );
        quantities_to_read.insert("psffit.amplitude");
        quantities_to_read.insert("psffit.magnitude_1adu");
        quantities_to_read.insert("psffit.mag");
        quantities_to_read.insert("projsrc.x");
        quantities_to_read.insert("projsrc.y");
        std::string psf_model=psf_data.get<std::string>(
            "psffit.model",
            "",
            IO::translate_string
        );
        if(psf_model=="") throw Error::CommandLine(
                "Input PSF map file does not define a PSF model."
        );
        if(psf_model=="bicubic" || psf_model=="zero")
            return new PSF::PiecewiseBicubicMap(psf_data, max_aperture + 1.0);
        else {
            assert(psf_model=="sdk");
            return new PSF::EllipticalGaussianMap(psf_data);
        }
    }

    ///List all the keys in the given data tree to stdout.
    void print_data_tree_keys(
        const boost::property_tree::basic_ptree<
            std::basic_string<char>,
            boost::any
        > &data_tree,
        const std::string &indent = ""
    )
    {
        for(
            IO::H5IODataTree::const_iterator it = data_tree.begin();
            it != data_tree.end();
            ++it
        ) {
            std::cout << indent << it->first << std::endl;
            print_data_tree_keys(it->second, indent+"\t");
        }
    }

} //End SubPixPhot namespace.

