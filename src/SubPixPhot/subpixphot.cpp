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
#include "../IO/SubPixHDF5File.h"
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
#include <H5Cpp.h>
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
        IO::SubPixHDF5File psfmap_file(psfmap_filename.c_str(),
                                       H5F_ACC_RDONLY);
        psfmap_file.read(quantities_to_read.begin(),
                         quantities_to_read.end(),
                         psf_data,
                         false);
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

    ///\brief Reads the positions and background (if available) of the source to
    ///photometer.
    ///
    ///If backgrounds are not available in the input file they are measured.
    void read_input_sources(
        ///The name of the HDF5 file to read the sources from.
        const std::string &source_filename,

        ///The data tree to fill.
        IO::H5IODataTree &data_tree)
    {
        std::set<std::string> quantities_to_read;
        quantities_to_read.insert("projsrc.x");
        quantities_to_read.insert("projsrc.y");
        quantities_to_read.insert("bg.value");
        quantities_to_read.insert("bg.error");
        quantities_to_read.insert("psffit.flux");
        quantities_to_read.insert("psffit.amplitude");
        quantities_to_read.insert("psffit.magnitude_1adu");
        quantities_to_read.insert("psffit.mag");
        quantities_to_read.insert("psffit.model");
        IO::SubPixHDF5File source_list_file(source_filename.c_str(),
                                            H5F_ACC_RDONLY);
        source_list_file.read(quantities_to_read.begin(),
                              quantities_to_read.end(),
                              data_tree,
                              false);
        PSF::fill_psf_amplitudes(data_tree);
    }

} //End SubPixPhot namespace.

///Perform the photometry specified through the command line.
int main(int argc, char **argv)
{
#ifndef DEBUG
//	try {
#endif
    SubPixPhot::Config options(argc, argv);
    if(!options.proceed()) return 1;
#ifdef TRACK_PROGRESS
    std::cerr << "Parsed command line." << std::endl;
#endif
    IO::FitsImage<double> image(options["io.image"].as<std::string>());
    IO::FitsImage<double> subpix_map(options["io.subpix"].as<std::string>());
    Core::SubPixelMap default_subpix_map(1, 1, "constant");
    default_subpix_map(0, 0)=1;
    if(options["io.sources"].as<std::string>()=="-")
        throw Error::CommandLine("Reading source from stdin not "
                                 "supported.");

    Core::RealList apertures=options["ap.aperture"].as<Core::RealList>();
    apertures.sort();

    IO::H5IODataTree data_tree(argc, argv, SUB_PIX_PHOT_VERSION, options);

    PSF::Map *psf_map=SubPixPhot::get_psf_map(
        options["io.psfmap"].as<std::string>(),
        apertures.back(),
        data_tree
    );


    SubPixPhot::read_input_sources(options["io.sources"].as<std::string>(),
                                   data_tree);
    bool output_background=false;
    if(!data_tree.get_optional<boost::any>("bg.value")) {
        SubPixPhot::add_background_measurements(
            options["bg.value"].as<Background::Annulus>(),
            image,
            data_tree
        );
        output_background=true;
    }

    if(options["io.subpix"].as<std::string>()=="") {
        Core::SubPixelCorrectedFlux<Core::SubPixelMap> measure_flux(
            image,
            default_subpix_map,
            options["ap.const-error"].as<double>(),
            apertures,
            options["gain"].as<double>()
        );
        SubPixPhot::add_flux_measurements(
            *psf_map,
            measure_flux,
            options["magnitude-1adu"].as<double>(),
            data_tree
        );
    } else {
        Core::SubPixelCorrectedFlux< IO::FitsImage<double> >
            measure_flux(
                image,
                subpix_map,
                options["ap.const-error"].as<double>(),
                apertures,
                options["gain"].as<double>()
            );
        SubPixPhot::add_flux_measurements(
            *psf_map,
            measure_flux,
            options["magnitude-1adu"].as<double>(),
            data_tree
        );
    }

    if(!output_background) {
        data_tree.erase("bg");
    }
    data_tree.erase("psffit");
    data_tree.erase("projsrc");

    IO::SubPixHDF5File *file;
    try {
        file=new IO::SubPixHDF5File(
            options["io.output"].as<std::string>().c_str(),
            H5F_ACC_RDWR
        );
    } catch(H5::FileIException) {
        try {
            file=new IO::SubPixHDF5File(
                options["io.output"].as<std::string>().c_str(),
                H5F_ACC_TRUNC
            );
        } catch(H5::FileIException &ex) {
            ex.printErrorStack();
            return 1;
        }
    }
    file->write(data_tree, false);
    file->close();
    delete file;
    delete psf_map;
    return 0;
#ifndef DEBUG
    /*	} catch(Error::General &ex) {
        std::cerr << ex.what() << ": " << ex.get_message() << std::endl;
        return 2;
        }*/
#endif
}
