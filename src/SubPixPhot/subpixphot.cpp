/** \file SubPixPhot.h
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
#include <H5Cpp.h>
#include <string>
#include <list>
#include <ctime>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>


const std::string SUB_PIX_PHOT_VERSION="$Revision: $";

namespace SubPixPhot {

    void Config::describe_options()
    {
        _hidden.add_options()
            (
                "io.image", 
                opt::value<std::string>(),
                "Reduced file to fit the PSF of."
            )
            (
                "io.sources,l",
                opt::value<std::string>(),
                "HDF5 file listing the sources to perform photometry for. "
                "Typically produced by FitPSF."
            );
        _positional.add("io.image", 1);
        _positional.add("io.sources", 1);
        opt::options_description io_options("Input/Output options");
        io_options.add_options()
            (
                "io.psfmap",
                opt::value<std::string>(),
                "HDF5 file defining the PSF map. If not specified, defaults "
                "to the input sources file."
            )
            (
                "io.output,o",
                opt::value<std::string>(),
                "File to place the output in HDF5 format. If not specified, "
                "defaults to the input sources file."
            )
            (
                "io.subpix,s",
                opt::value<std::string>()->default_value(""),
                "FITS file representing the sensitivity structure of a "
                "pixel. If not specified, uniform sensitivity is assumed."
            )
            (
                "io.output-quantities",
                opt::value<Core::ColumnList>(),
                "A comma separated list of quantities to output. Recognized "
                "values: id, x, y, S, D, K, A|amp, flux, flux_err, mag, "
                "mag_err, flag, bg, bg_err, enabled, sn, npix, nbgpix"
#ifdef DEBUG
                ", time"
#endif
                ". In text output a group of consecutive columns that depend"
                " on aperture are repeated automatically."
            );
        opt::options_description apphot_options(
            "Options configuring how aperture photometry is performed.");
        apphot_options.add_options()
            (
                "ap.aperture",
                opt::value<Core::RealList>(),
                "Comma separated list of apertures to use."
            )
            (
                "ap.const-error",
                opt::value<double>(),
                "A value to add to the error estimate of a pixel (intended "
                "to represent things like readout noise, truncation noise "
                "etc.). "
            );
        opt::options_description generic_options(
            "Options relevant for more than one component.");
        generic_options.add_options()
            (
                "gain,g",
                opt::value<double>(),
                "How many electrons are converted to 1ADU."
            )
            (
                "magnitude-1adu",
                opt::value<double>(),
                "The magnitude that corresponds to a flux of 1ADU on the "
                "input image."
            );
        opt::options_description background_options(
            "Options defining how to measure the background behing the "
            "sources.");
        background_options.add_options()
            (
                "bg.annulus,b",
                opt::value<Background::Annulus>(),
                "Specifies that an annulus with the given inner radius "
                "centered around the source should be used to estimate the "
                "background and its error."
            )
            (
                "bg.min-pix",
                opt::value<unsigned>(),
                "If a source's background is based on less than this many "
                "pixels, the source is excluded from the fit."
            );
        _cmdline_config.add(io_options)
            .add(generic_options)
            .add(apphot_options)
            .add(background_options);
    }

    void SubPixPhot::Config::apply_fallbacks()
    {
        std::stringstream fallbacks;
        fallbacks 
            << "[io]" << std::endl
            << "psfmap = " << (*this)["io.sources"].as<std::string>()
            << std::endl
            << "output = " << (*this)["io.sources"].as<std::string>()
            << std::endl;
        fallbacks.seekg(0);
        opt::store(
            opt::parse_config_file(fallbacks,
                                   _cmdline_config),
            *this
        );
        opt::notify(*this);
    }

    void SubPixPhot::Config::check_consistency()
    {
        if((*this)["ap.aperture"].as<Core::RealList>().size() == 0)
            throw Error::CommandLine("No apertures defined!");
        for(
                Core::RealList::const_iterator
                    i = (*this)["ap.aperture"].as<Core::RealList>().begin();
                i != (*this)["ap.aperture"].as<Core::RealList>().end();
                ++i
        )
            if(*i < 0) 
                throw Error::CommandLine("Negative aperture encountered!");
        if((*this)["gain"].as<double>()<0)
            throw Error::CommandLine("Negative gain specified!");
        assert(count("io.psfmap"));
        assert(count("io.output"));
    }

    std::string SubPixPhot::Config::usage_help(
        const std::string &prog_name
    ) const
    {
        return "Usage: "
               +
               prog_name
               +
               " [OPTION ...] <FITS image> <HDF5 sources>";
    }

    ///Measures the fluxes of the sources using aperture photometry.
    template<class FLUX_MEASURER>
    void add_flux_measurements(
            ///The SPF map to use.
            const PSF::Map &psf_map, 

            ///The object that will measure the fluxes.
            FLUX_MEASURER &measure_flux,

            ///The magnitude that corresponds to a flux of 1ADU
            double mag_1adu,

            ///The data tree to add the fluxes to. It must contain the source
            ///positions and backgrounds on input.
            IO::H5IODataTree &data_tree)
    {
        unsigned num_apertures = measure_flux.number_apertures();

        PSF::MapSourceContainer psfmap_sources(data_tree, num_apertures);
        unsigned num_sources = psfmap_sources.size();

        IO::OutputArray<double>
            background(data_tree.get<boost::any>("bg.value")),
            background_error(data_tree.get<boost::any>("bg.error"));

        std::vector< std::vector<double>* > magnitudes(num_apertures),
                                            magnitude_errors(num_apertures);
        std::vector< std::vector<unsigned>* > flags(num_apertures);
        for(unsigned ap_index = 0; ap_index < num_apertures; ++ap_index) {
            magnitudes[ap_index] = new std::vector<double>(num_sources);
            magnitude_errors[ap_index] = 
                new std::vector<double>(num_sources);
            flags[ap_index] = new std::vector<unsigned>(num_sources);
        }
                                           
        PSF::MapSourceContainer::const_iterator 
            psfmap_src_iter = psfmap_sources.begin();

        unsigned x_resolution = measure_flux.image().x_resolution(),
                 y_resolution = measure_flux.image().y_resolution();
        for(
                size_t source_index = 0;
                source_index < num_sources;
                ++source_index
        ) {
#ifdef TRACK_PROGRESS
            clock_t t1_clock = std::clock();
            time_t t1_time = std::time(0);
#endif
            std::valarray<Core::Flux> measured_flux_values(num_apertures);
            double x = psfmap_src_iter->x(),
                   y = psfmap_src_iter->y();
            if(
                    x < 0 || x > x_resolution || y < 0 || y > y_resolution
                    ||
                    std::isnan(psfmap_src_iter->background().value())
            )
                measured_flux_values.resize(
                    num_apertures,
                    Core::Flux(Core::NaN, Core::NaN, Core::BAD)
                );
            else {
                PSF::PSF *psf = psf_map(
                    psfmap_src_iter->expansion_terms(),
                    psfmap_src_iter->background().value()
                );
                measured_flux_values = measure_flux(
                        x,
                        y,
                        *psf,
                        background[source_index],
                        background_error[source_index]
                );
                delete psf;
            }
            for(unsigned ap_index = 0; ap_index < num_apertures; ++ap_index){
                const Core::Flux &measured = measured_flux_values[ap_index];
                (*(magnitudes[ap_index]))[source_index] = magnitude(
                        measured.value(),
                        mag_1adu
                );
                (*(magnitude_errors[ap_index]))[source_index] = 
                    magnitude_error(
                        measured.value(),
                        measured.error()
                    );
                (*(flags[ap_index]))[source_index] = measured.flag();
            }
#ifdef TRACK_PROGRESS
            clock_t t2_clock = std::clock();
            time_t t2_time = std::time(0);
            std::cerr << "Source " << source_index << " took "
                      << (t2_clock == clock_t(-1)
                          ? std::difftime(t2_time, t1_time)
                          : double(t2_clock - t1_clock) / CLOCKS_PER_SEC)
                      << " seconds to photometer" << std::endl;
#endif
            ++psfmap_src_iter;
        }

        IO::TranslateToAny< std::vector<double> > double_trans;
        IO::TranslateToAny< std::vector<unsigned> > unsigned_trans;
        for(unsigned ap_index = 0; ap_index < num_apertures; ++ap_index) {
            std::ostringstream key;
            key << "apphot.mag." << ap_index;
            data_tree.put(key.str(), magnitudes[ap_index], double_trans);
            key.str("");
            key << "apphot.mag_err." << ap_index;
            data_tree.put(key.str(),
                          magnitude_errors[ap_index],
                          double_trans);
            key.str("");
            key << "apphot.quality." << ap_index;
            data_tree.put(key.str(), flags[ap_index], unsigned_trans);
        }
    }

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
