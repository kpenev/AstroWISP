/**\file
 *
 * \brief The definitions of the functions declared in DataTreeCalculations.h
 *
 * \ingroup IO
 */

#include "DataTreeCalculations.h"

namespace PSF {

    void fill_psf_fluxes(IO::H5IODataTree &data,
                         const std::string &data_tree_image_id)
    {
        std::string tree_suffix = (data_tree_image_id.empty()
                                   ? ""
                                   : "." + data_tree_image_id);
        if(data.get_optional<boost::any>("psffit.flux" + tree_suffix))
            return;
        double mag_1adu=data.get<double>("psffit.magnitude_1adu",
                                         Core::NaN,
                                         IO::translate_double);
        IO::OutputArray<double>
            mag_array(data.get<boost::any>("psffit.mag" + tree_suffix));
        Eigen::Map< const Eigen::ArrayXd > mag_eigen(mag_array.data(),
                                                      mag_array.size(),
                                                      1);
        Eigen::ArrayXd *flux=new Eigen::ArrayXd(Eigen::exp(
                std::log(10.0)*2.0/5.0*(mag_1adu-mag_eigen)
        ));
        data.put("psffit.flux" + tree_suffix,
                 flux,
                 IO::TranslateToAny< Eigen::ArrayXd >());
    }

    void fill_psf_amplitudes(IO::H5IODataTree &data,
                             const std::string &data_tree_image_id)
    {
        std::string tree_suffix = (data_tree_image_id.empty()
                                   ? ""
                                   : "." + data_tree_image_id);
        typedef Eigen::Map< const Eigen::ArrayXd > MapType;
        if(data.get_optional<boost::any>("psffit.amplitude" + tree_suffix))
            return;
        fill_psf_fluxes(data, data_tree_image_id);
        std::string psf_model=data.get<std::string>("psffit.model",
                                                    "",
                                                    IO::translate_string);
        if(psf_model=="bicubic" || psf_model=="zero") {
            data.put("psffit.amplitude" + tree_suffix,
                     data.get<boost::any>("psffit.flux"));
            return;
        }
        assert(psf_model=="sdk");
        IO::OutputArray<double> fluxes_array(
            data.get<boost::any>("psffit.flux" + tree_suffix)
        );
        size_t num_sources=fluxes_array.size();
        MapType flux(fluxes_array.data(), num_sources, 1),
                s(NULL, num_sources, 1),
                d(NULL, num_sources, 1),
                k(NULL, num_sources, 1);

        MapSourceContainer psfmap_sources(data, 1);
        assert(psfmap_sources.size() == num_sources);
        if(psfmap_sources.front().expansion_terms().size() == 0) {
            IO::OutputArray<double>
                s_array(data.get<boost::any>("psffit.s")),
                d_array(data.get<boost::any>("psffit.d")),
                k_array(data.get<boost::any>("psffit.k"));
            assert(s_array.size()==num_sources);
            assert(d_array.size()==num_sources);
            assert(k_array.size()==num_sources);
            new (&s) MapType(s_array.data(), num_sources, 1);
            new (&d) MapType(d_array.data(), num_sources, 1);
            new (&k) MapType(k_array.data(), num_sources, 1);
        } else {

            EllipticalGaussianMap psf_map(data);
            Eigen::ArrayXd s_array(num_sources),
                           d_array(num_sources),
                           k_array(num_sources);
            MapSourceContainer::const_iterator
                src_i = psfmap_sources.begin();
            for(size_t i=0; i<num_sources; ++i) {
                EllipticalGaussian *psf = psf_map.get_psf(
                    src_i->expansion_terms()
                );
                s_array[i]=psf->s();
                d_array[i]=psf->d();
                k_array[i]=psf->k();
                ++src_i;
            }
            new (&s) MapType(s_array.data(), num_sources, 1);
            new (&d) MapType(d_array.data(), num_sources, 1);
            new (&k) MapType(k_array.data(), num_sources, 1);
        }
        Eigen::ArrayXd *amplitude=new Eigen::ArrayXd(
                flux*Eigen::sqrt(s.square()-d.square()-k.square())/(2.0*M_PI)
        );
        data.put("psffit.amplitude" + tree_suffix,
                 amplitude,
                 IO::TranslateToAny< Eigen::ArrayXd >());
    }

} //End IO namespace.
