/**\file
 * \brief The implementation of some of the Variance methods and the methods
 * converting between a sub-pixel sensitivity map and a unit cube.
 *
 * \ingroup FitSubpix
 */
#include "FittingUtil.h"

void Variance::fill_image_w(const FitsImage &frame,
                            const H5IODataTree &source_tree,
                            SubPixelCorrectedFlux<SubPixelMap> &flux)
{
    OutputArray<double> x(source_tree.get<boost::any>("projsrc.x")),
                        y(source_tree.get<boost::any>("projsrc.y")),
                        amp(source_tree.get<boost::any>("psffit.amplitude")),
                        bg(source_tree.get<boost::any>("bg.value")),
                        bg_err(source_tree.get<boost::any>("bg.error")),

    PSFMap *psf_map;
	if(psf_model=="") throw Error::CommandLine(
			"Input PSF map file does not define a PSF model."
	);
	if(psf_model=="bicubic") 
            psf_map=new PiecewiseBicubicPSFMap(source_tree, max_aperture);
	else {
		assert(psf_model=="sdk");
		psf_map=new EllipticalGaussianPSFMap(source_tree);
	}
    hsize_t nsources=x.size();
    assert(y.size()==nsources);
    assert(bg.size()==nsources);
    assert(bg_err.size()==nsources);

    for(size_t s=0; s<nsources; ++s) {
        PSF *psf=psf_map(x[s], y[s], bg[s]/amp[s]);
        flux(x[s], y[s], *psf, bg[s], bg_err[s], true);
        __w[s][t]=flux.get_subpix_scaling(0);
    }
}

void Variance::calculate_wb(const StringList &frame_fnames,
                            const StringList &source_fnames,
                            double aperture)
{
	assert(frame_fnames.size()==source_fnames.size());
	std::list<double> ap_list(1, aperture);
	SubPixelMap dummy_map(__num_x_subpix, __num_y_subpix, "dummy");
	SubPixelCorrectedFlux<SubPixelMap> flux(dummy_map, 0, ap_list);
    hsize_t nsources = 0;
    size_t frame_count = frame_fnames.size();
    StringList::const_iterator frame_fname_i=frame_fnames.begin(),
                               source_fname_i=source_fnames.begin();
    OutputArray<std::string> first_names;
    OutputArray<unsigned first_fields, first_sources;

    std::ostringstream source_mismatch_msg;
    source_mismatch_msg << "Different source lists found in '"
                        << frame_fnames[0]
                        << "' and '";

	for(size_t t=0; t<frame_count; t++) {
		FitsImage<double> frame(*frame_fname_i);
		flux.set_image(frame);
        H5IODataTree source_tree;
        read_psf_file(*source_fname_i, source_tree);
        if(t) {
            OutputArray<std::string> 
                these_names(source_tree.get<boost:any>("projsrc.srcid.name"));
            OutputArray<unsigned>
                these_fields(
                        source_tree.get<boost::any>("projsrc.srcid.field")
                ),
                these_sources(
                        source_tree.get<boost::any>("projsrc.srcid.source")
                );
            if(
                    these_sources!=first_sources 
                    || these_fields!=first_fields 
                    || these_names!=first_names
            ) {
                source_mismatch_msg << *frame_fname_i << "'!";
                throw Error::IO(msg.str());
            }
        } else {
            first_fields.parse(
                    source_tree.get<boost::any>("projsrc.srcid.field")
            );
            first_sources.parse(
                    source_tree.get<boost::any>("projsrc.srcid.source")
            );
            first_names.parse(
                    source_tree.get<boost::any>("projsrc.srcid.name")
            );
            assert(first_sources.size()==first_fields.size());
            nsources=first_names.size();
            if(nsources==0) nsources=first_sources.size();

			__w.resize(nsources);
			__ws_inverse.resize(nsources);
            __xi.resize(nsources, frame_count);
            __f.resize(nsources, frame_count);
            for(hsize_t s=0; s<nsources; ++s) {
				__w[s].resize(frame_count);
				__ws_inverse[s].resize(frame_count);
            }
		}
        fill_image_w(frame, source_tree, flux);

        ++frame_fname_i;
        ++source_fname_i;
	}
}

void Variance::calculate_lambda(std::valarray<Eigen::ArrayXXd> &lambda, 
		std::valarray<Eigen::ArrayXd> &lambda_bar,
		unsigned num_subpix)
{
	unsigned max_s=__w.size(), max_t=__w[0].size();
	for(unsigned k=0; k<num_subpix; k++) {
		lambda[k].resize(max_s, max_t);
		for(unsigned s=0; s<max_s; s++)
			for(unsigned t=0; t<max_t; t++)	
				lambda[k](s,t)=__w[s][t].col(k).dot(
					__ws_inverse[s][t].array().square().matrix());
		lambda_bar[k]=(__xi.array()*lambda[k]).rowwise().sum();
	}
}

void Variance::calculate_omega(
    Eigen::ArrayXXd&    omega,
    unsigned            k,
    unsigned            l
)
{
	unsigned max_s=__w.size(), max_t=__w[0].size();
	for(unsigned s=0; s<max_s; s++)
		for(unsigned t=0; t<max_t; t++) {
			omega(s,t)=2.0*(__w[s][t].col(k).array()*
					__w[s][t].col(l).array()*
					__ws_inverse[s][t].array().cube()).sum();
		}
}

void Variance::initialize(const StringList &frame_fnames,
                          const StringList &source_fnames,
                          double aperture,
                          unsigned subpix_x_res,
                          unsigned subpix_y_res
)
{
	__num_x_subpix=subpix_x_res; 
	__num_y_subpix=subpix_y_res;
	calculate_wb(frame_fnames, source_fnames, aperture);
	calculate_q();
	calculate_xi();
}

Eigen::MatrixXd calculate_jacobian(const Eigen::VectorXd &x)
{
	Eigen::MatrixXd jac;
	unsigned N=x.size();
	jac.setZero(N, N+1);
	double range=N+1;
	for(unsigned k=0; k<N+1; k++) {
		double sk=(k<N ? range*(1.0-std::pow(1.0-x(k), 1.0/(N-k))) :range);
		range-=sk;
		for(unsigned l=0; l<k; l++)
			jac(l,k)=-sk/((1.0-x(l))*(N-l));
		if(k<N) jac(k,k)=range/((1.0-x(k))*(N-k));
	}
	return jac;
}

std::valarray<Eigen::MatrixXd> calculate_d2s(const Eigen::VectorXd &x)
{
	unsigned N=x.size();
	std::valarray<Eigen::MatrixXd> result(Eigen::MatrixXd::Zero(N,N), N+1);
	double range=N-1;
	for(unsigned i=0; i<N+1; i++) {
		unsigned Nmi=N-i;
		double xi=(Nmi ? x(i) : 1), 
			   si=(Nmi ? range*(1.0-std::pow(1.0-xi, 1.0/(Nmi))) : range);
		range-=si;
		for(unsigned k=0; k<i; k++) {
			double xk=x(k);
			unsigned Nmk=N-k;
			for(unsigned l=0; l<k; l++)
				result[i](k,l)=result[i](l,k)=
					si/((N-l)*Nmk*(1.0-x(l))*(1.0-xk));
			result[i](k,k)=si/(Nmk*std::pow(1.0-xk, 2))*(1.0/Nmk-1.0);
			if(i<N) result[i](k,i)=result[i](i,k)=-range/(
					Nmi*Nmk*(1.0-xi)*(1.0-xk));
		}
		if(i<N) result[i](i,i)=range/(Nmi*std::pow(1.0-xi,2))*(1.0-1.0/Nmi);
	}
	return result;
}

void ChiSquared::read_psf_file(const std::string &filename,
                               unsigned image_index,
                               H5IODataTree &data_tree)
{
    std::set<std::string> quantities_to_read(
        PiecewiseBicubicPSFMap::required_data_tree_quantities()
	);
    quantities_to_read.insert(
        EllipticalGaussianPSFMap::required_data_tree_quantities().begin(),
        EllipticalGaussianPSFMap::required_data_tree_quantities().end()
	);
	quantities_to_read.insert("projsrc.x");
	quantities_to_read.insert("projsrc.y");
	quantities_to_read.insert("projsrc.srcid.field");
	quantities_to_read.insert("projsrc.srcid.source");
	quantities_to_read.insert("projsrc.srcid.name");
	quantities_to_read.insert("psffit.amplitude");
	quantities_to_read.insert("bg.value");
	quantities_to_read.insert("bg.error");
	SubPixHDF5File psf_file(filename.c_str(), H5F_ACC_RDONLY);
	psfmap_file.read(quantities_to_read.begin(),
					 quantities_to_read.end(),
					 data_tree,
					 false);
}
