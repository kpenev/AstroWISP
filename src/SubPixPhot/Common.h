/**\file
 *
 * \brief Declare functions used by both the C-interface and the command line
 * interface.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SUBPIX_PHOT_COMMON_H
#define __SUBPIX_PHOT_COMMON_H

#include "../IO/OutputArray.h"
#include "../IO/H5IODataTree.h"
#include "../IO/TranslateToAny.h"
#include "../PSF/Map.h"
#include "../PSF/MapSourceContainer.h"
#include "../PSF/PSF.h"
#include "../Core/Flux.h"

namespace SubPixPhot {

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
        IO::H5IODataTree &data_tree,

        ///See same name argument to
        ///PSF::MapSourceContainer::MapSourceContainer()
        const std::string &data_tree_image_id=""
    )
    {
        std::string tree_suffix = (data_tree_image_id.empty()
                                   ? ""
                                   : "." + data_tree_image_id);
        std::cerr << "tree suffix = " << tree_suffix << std::endl;
        unsigned num_apertures = measure_flux.number_apertures();

        PSF::MapSourceContainer psfmap_sources(data_tree,
                                               num_apertures,
                                               data_tree_image_id);
        unsigned num_sources = psfmap_sources.size();

        std::cerr << "Reading background" << std::endl;
        IO::OutputArray<double>
            background(
                data_tree.get<boost::any>("bg.value" + tree_suffix)
            ),
            background_error(
                data_tree.get<boost::any>("bg.error" + tree_suffix)
            );

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
#ifndef NDEBUG
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
            ) {
                measured_flux_values.resize(
                    num_apertures,
                    Core::Flux(Core::NaN, Core::NaN, Core::BAD)
                );
            } else {
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
#ifndef NDEBUG
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
            key << "apphot.mag" << tree_suffix << "." << ap_index;
            data_tree.put(key.str(), magnitudes[ap_index], double_trans);
            key.str("");
            key << "apphot.mag_err" << tree_suffix << "." << ap_index;
            data_tree.put(key.str(),
                          magnitude_errors[ap_index],
                          double_trans);
            key.str("");
            key << "apphot.quality" << tree_suffix << "." << ap_index;
            data_tree.put(key.str(), flags[ap_index], unsigned_trans);
        }
    }

} //End SubPixPhot namespace.

#endif
