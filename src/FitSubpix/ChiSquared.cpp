/// \ingroup FitSubpix 

/**\file
 * \brief The implementation of the non-trivial ChiSquared class methods.
 */

#include "ChiSquared.h"
#include <iostream>
#include <iomanip>

void ChiSquared::verify_sources(const H5IODataTree &this_tree
                                const std::string &psf_fname)
{
    OutputArray<std::string>
        these_names(this_tree.get<boost::any>("projsrc.srcid.name"));

    std::ostringstream msg;
    msg << "Different source list found in '"
        << psf_fname 
        << "' than in the first source list file!";

    if(these_names.size()==0) {
        OutputArray<unsigned> 
            these_fields(this_tree.get<boost::any>("projsrc.srcid.field")),
            these_sources(this_tree.get<boost::any>("projsrc.srcid.source")),
            first_fields(
                    __data_trees[0].get<boost::any>("projsrc.srcid.field")
            ),
            first_sources(
                    __data_trees[0].get<boost::any>("projsrc.srcid.source")
            );

        if(these_sources!=first_sources || these_fields!=first_fields)
            throw Error::IO(msg.str());
    } else {
        OutputArray<std::string> first_names(
                __data_trees[0].get<boost::any>("projsrc.srcid.name")
        );
        if(these_names!=first_names) throw Error::IO(msg.str());
    }
}

void ChiSquared::measure_image_fluxes(unsigned image_index,
                                      const SubPixelMap &subpix_map,
                                      const std::list<double> &aperture_list)
{
    SubPixelCorrectedFlux<SubPixelMap> measure_flux(__images[image_index],
                                                    subpix_map,
                                                    0.0,
                                                    aperture_list);
    OutputArray<double> 
        x(__data_trees[image_index].get<boost::any>("projsrc.x")),
        y(__data_trees[image_index].get<boost::any>("projsrc.y")),
        amplitude(
                __data_trees[image_index].get<boost::any>("psffit.amplitude")
        ),
        background(__data_trees[image_index].get<boost::any>("bg.value")),
        background_err(
                __data_trees[image_index].get<boost::any>("bg.error")
        );
    hsize_t num_sources=x.size();
    assert(y.size()==num_sources);
    assert(background.size()==num_sources);
    assert(background_err.size()==num_sources);

	std::string psf_model=__data_trees[image_index].get<std::string>(
        "psffit.model", "", translate_string
    );
    PSFMap *psf_map;
	if(psf_model=="") throw Error::CommandLine(
			"Input PSF map file does not define a PSF model."
	);
	if(psf_model=="bicubic") 
        psf_map=new PiecewiseBicubicPSFMap(__data_trees[image_index],
                                           max_aperture);
	else {
		assert(psf_model=="sdk");
		psf_map=new EllipticalGaussianPSFMap(__data_trees[image_index]);
	}

	double bg_scaling=M_PI*std::pow(aperture, 2);
    unsigned flux_index=0;
    for(unsigned src_ind=0; src_ind<num_sources; ++src_ind) {
        if(!commons_source[src_ind]) continue;
        PSF *psf=psf_map(x[src_ind],
                         y[src_ind],
                         background[src_ind]/amplitude[src_ind]);
        Flux flux=measure_flux(x[src_ind], y[src_ind],  *psf, true)[0];
        __flux_values[flux_index][image_index]=
            flux.value() - background[src_ind]*bg_scaling;
        if(flux.flag()==GOOD) 
            __flux_errors[flux_index][image_index]=std::sqrt(
                    std::pow(flux.error(), 2)
                    +
                    std::pow(background_err[src_ind]*bg_scaling, 2)
            );
        else __flux_errors[flux_index][image_index]=
                std::numeric_limits<double>::infinity();
        delete psf;
        ++flux_index;
    }
    delete psf_map;
}

double ChiSquared::chi2() const
{
	std::valarray<double> averages(__flux_values.size());
	std::valarray<double> individual_chi2(__flux_values.size());
	std::valarray<double> expected_var(__flux_values.size());
	for(size_t i=0; i<__flux_values.size(); i++) {
//		std::valarray<double> weights(1.0, __flux_errors[i].size());//1.0/__flux_errors[i];
//		weights/=weights.sum();
//		averages[i]=(__flux_values[i]*weights).sum();
		averages[i]=__flux_values[i].sum()/__flux_values[i].size();
		individual_chi2[i]=std::pow((__flux_values[i]-averages[i]), 2).sum()/__flux_values[i].size();
//		individual_chi2[i]=std::pow((__flux_values[i]-averages[i])*weights, 2).sum();
//		expected_var[i]=std::sqrt(std::pow(__flux_errors[i]*weights,2).sum());
	}
	return individual_chi2.sum()/__flux_values.size();
	return (individual_chi2/expected_var).sum();
}

void ChiSquared::initialize(
				const StringList &frame_filenames,
				const StringList &psf_filenames,
				double aperture)
{
    unsigned num_images=frame_filenames.size();
	assert(psf_filenames.size()==num_images);

	__images.assign(frame_filenames.begin(), frame_filenames.end());
	__aperture=aperture;
	__max_exp_coef=max_exp_coef;
    __data_trees.resize(num_images);

	double bg_scaling=M_PI*std::pow(aperture, 2);
	for(size_t i=0; i<frame_filenames.size(); i++) {
        read_psf_file(psf_filenames[i], __data_trees[i]);
		if(i) verify_sources(__data_trees[i], psf_filenames[i]);
	}
	__flux_values.resize(__sources[0].size());
	__flux_errors.resize(__sources[0].size());
	for(size_t i=0; i<__flux_values.size(); i++) {
		__flux_values[i].resize(frame_filenames.size());
		__flux_errors[i].resize(frame_filenames.size());
	}
}

double ChiSquared::operator()(const SubPixelMap &subpix_map) const
{
	std::list<double> aperture_list(1, __aperture);
	SubPixelMap normalized_map(subpix_map.x_resolution(), 
					           subpix_map.y_resolution(),
                               "copy");
	double sum=0.0;
	for(unsigned long y=0; y<subpix_map.y_resolution(); y++)
		for(unsigned long x=0; x<subpix_map.x_resolution(); x++)
			sum+=subpix_map(x,y);
	sum/=subpix_map.x_resolution()*subpix_map.y_resolution();
	for(unsigned long y=0; y<subpix_map.y_resolution(); y++)
		for(unsigned long x=0; x<subpix_map.x_resolution(); x++)
			normalized_map(x,y)=subpix_map(x,y)/sum;
	for(size_t image_index=0; image_index<__images.size(); image_index++)
        measure_image_fluxes(image_index, normalized_map, aperture_list);
	return chi2();
}
