#ifndef __FITTING_UTIL_H
#define __FITTING_UTIL_H

/** \file
 * \brief Various utilities useful for more than one fitting method of
 * FitSubpix. 
 *
 * \ingroup FitSubpix
 */

#include "SubPixelCorrectedFlux.h"
#include "SubPixelMap.h"
#include "PhotColumns.h"
#include "Background.h"
#include "CommandLineUtil.h"
#include "SubPixPhotIO.h"
#include "Eigen/Dense"
#include <valarray>
#include <cassert>
#include <vector>
#include <string>

///\brief The Variance and its first two derivatives with respect to subpixel
///sensitivities.
///
///Allows for arbitrary weights of the various stars and observations.
///Details of how everything is calculated are given
/// <a href="NRFitting.pdf">here</a> The variable names are the same as used
///there.
///
///\ingroup FitSubpix
class LIB_LOCAL Variance {
private:
	///A vector of the relative weights of the various images.
	Eigen::RowVectorXd __q;
	
	///Rows are the relative weights of different observations of a star.
	Eigen::MatrixXd __xi;

	///\brief The sub-pixel independent scaling matrices used to compute
	///fluxes for a given subpixel map.
	std::valarray< std::valarray<Eigen::MatrixXd> > __w;

	///\brief Each entry is the inverse of the corresponding entry in the
	/// \f$ w \cdot s \f$ vector for each observation of each star.
	std::valarray< std::valarray<Eigen::VectorXd> > __ws_inverse;

	///\brief Rows are the fluxes of a source on all frames corrected with
	///the sub-pixel structure.
	Eigen::MatrixXd __f;
	
	unsigned __num_x_subpix, ///< The number of x splittings of a pixel
			 __num_y_subpix; ///< The number of y splittings of a pixel

	///The minimum subdivision size allowed when integrating PSFs
	double __max_exp_coef;

    ///Fills a slice of __w corresponding to a single image.
    void fill_image_w(
            ///The image to process.
            const FitsImage<double> &frame,

            ///The psf fitted sources in the frame.
            const H5IODataTree &source_tree,

            ///Object to measure aperture photometry flux from the image.
            ///Modified, because the image on which it will operate is set.
            SubPixelCorrectedFlux<SubPixelMap> &flux);

	///\brief Calculates the \f$w\f$ and \f$b\f$ matrices and the
	/// \f$\bar{b}\f$ vector.
	void calculate_wb(const StringList &frame_fnames,
			          const StringList &source_fnames,
                      double aperture);

	///\brief Calculates the weight of each image (\f$q\f$).
	///
	///The \f$w\f$ matrices must already be initialized.
	void calculate_q() {__q.setConstant(__w.size(), 1.0/__w.size());}

	///\brief Calculates the weight of each source in each image (\f$\xi\f$).
	///
	///The \f$w\f$ matrices must already be initialized.
	void calculate_xi()
	{__xi.setConstant(__w.size(), __w[0].size(), 1.0/__w[0].size());}

	///Fills #__ws_inverse and #__f with the appropriate values.
	template<typename EIGEN_TP>
		void calculate_ws_inverse_f(const Eigen::MatrixBase<EIGEN_TP> &s);


	///\brief Calculates all \f$\Lambda\f$ matrices.
	///
	///The array should have the correct size but the matrices are properly
	///resized.
	void calculate_lambda(std::valarray<Eigen::ArrayXXd> &lambda,
			std::valarray<Eigen::ArrayXd> &lambda_bar,
			unsigned num_subpix);

	///Calculates the \f$\Omega_{kl}\f$ matrix.
    void calculate_omega(
        Eigen::ArrayXXd&    omega,
        unsigned            k,
        unsigned            l
    );

public:
	///Default constructor.
	Variance() : __max_exp_coef(1) {};

	///\brief Initializes all members which do not depend on the sub-pixel
	///structure.
    Variance(
        ///The filenames of the frames used for fitting
        const StringList &frame_fnames,

        ///The filenames of the lists of sources in each frame.
        const StringList &source_fnames,

        ///The aperture for measuring fluxes.
        double aperture,

        ///The x resolution of the sub-pixel sensitivity map.
        unsigned subpix_x_res,

        ///The y resolution of the sub-pixel sensitivity map.
        unsigned subpix_y_res
    ) : __max_exp_coef(1) 
    {
        initialize(
            frame_fnames,
            source_fnames,
            aperture,
            subpix_x_res,
            subpix_y_res
        );
    }

	///\brief Initializes all members which do not depend on the sub-pixel
	///structure.
	///
	///See the documentation of the constructor (Variance()) for a
	///description of the parameters.
    void initialize(
        const StringList&   frame_fnames,
        const StringList&   source_fnames,
        double              aperture,
        unsigned            subpix_x_res,
        unsigned            subpix_y_res
    );

	///\brief Calculates the flux variance for a given subpixel sensitivity
	///vector.
	template<typename EIGEN_TP>
		double variance(
				///The subpixel sensitivity vector to use.
				const Eigen::MatrixBase<EIGEN_TP> &s)
		{double result; calculate(s, &result); return result;}

	///\brief Calculates the gradient of the variance for a given subpixel
	///sensitivity.
	template<typename EIGEN_TP>
		Eigen::VectorXd grad_variance(
				///The subpixel sensitivity vector to use.
				const Eigen::MatrixBase<EIGEN_TP> &s)
		{Eigen::VectorXd result; calculate(s, NULL, &result); return result;}

	///\brief Calculates the second derivative matrix of the variance for a
	///given subpixel sinsitivity.
	template<typename EIGEN_TP>
		Eigen::MatrixXd d2variance(
				///The subpixel sensitivity vector to use.
				const Eigen::MatrixBase<EIGEN_TP> &s)
		{
			Eigen::MatrixXd result; 
			calculate(s, NULL, NULL, &result); 
			return result;
		}

	///\brief Calculates the variance, its gradient and its second derivative
	///matrix.
	///
	///If any of the output arguments are NULL, they are not
	///calculated.
	template <typename EIGEN_TP> void calculate(
			///The sib-pixel sensitivities
			const Eigen::MatrixBase<EIGEN_TP> &s,

			///The variance of the flux corresponding to the given
			///sensitivities, if NULL it is not calculated.
			double *variance=NULL,

			///The gradient of the variance with respect to the
			///sensitivities, if NULL it is not calculated.
			Eigen::VectorXd *grad=NULL,

			///The second derivative of the variance with respect to the
			///sensitivities, if NULL it is not calculated.
			Eigen::MatrixXd *d2=NULL);
};

///Convert a vector of unit cube values to the corresponding sub-pixel map.
template<typename EIGEN_DERIVED>
Eigen::VectorXd x_to_s(const Eigen::DenseBase<EIGEN_DERIVED> &x_vec);

///Convert a sub-pixel map to a vector of unit cube values.
template<typename EIGEN_DERIVED>
Eigen::VectorXd s_to_x(const Eigen::DenseBase<EIGEN_DERIVED> &s_vec);

///Calculates the Jocobian matrix for the given x vector.
Eigen::MatrixXd calculate_jacobian(const Eigen::VectorXd &x);

///\brief Calculates the second derivative matrices of each sub-pixel
///sensitivity with respect to x
std::valarray<Eigen::MatrixXd> calculate_d2s(const Eigen::VectorXd &x);

///\brief Read data from a psf file to a data_tree.
void ChiSquared::read_psf_file(
        ///The name of the file to read.
        const std::string &filename,

        ///The data tree to fill.
        H5IODataTree &data_tree);

template<typename EIGEN_DERIVED>
Eigen::VectorXd x_to_s(const Eigen::DenseBase<EIGEN_DERIVED> &x_vec)
{
	assert(x_vec.cols()==1);
	unsigned long N=x_vec.size();
	double range=N+1;
	Eigen::VectorXd s(N+1);
	for(unsigned i=0; i<N; i++) {
		double cube_val=x_vec(i); 
		s(i)=range*(1.0-std::pow(1.0-cube_val, 1.0/(N-i)));
		range-=s(i);
	}
	s(N)=range;
	return s;
}

template<typename EIGEN_DERIVED>
Eigen::VectorXd s_to_x(const Eigen::DenseBase<EIGEN_DERIVED> &s_vec)
{
	assert(s_vec.cols()==1);
	int N=s_vec.size()-1;
	double range=N+1;
	Eigen::VectorXd x(N);
	for(int i=0; i<N; i++) {
		x(i)=1.0-std::pow(1.0-s_vec(i)/range, N-i);
		range-=s_vec(i);
	}
	return x;
}

template<typename EIGEN_TP>         
void Variance::calculate_ws_inverse_f(const Eigen::MatrixBase<EIGEN_TP> &s)
{
	size_t max_s=__w.size(), max_t=__w[0].size();
	for(size_t s_ind=0; s_ind<max_s; s_ind++)
		for(size_t t_ind=0; t_ind<max_t; t_ind++) {
			__ws_inverse[s_ind][t_ind]=(__w[s_ind][t_ind]*s).array().inverse().matrix();
			__f(s_ind,t_ind)=__ws_inverse[s_ind][t_ind].sum();
		}
}

template <typename EIGEN_TP> void Variance::calculate(
		const Eigen::MatrixBase<EIGEN_TP> &s, double *variance,
		Eigen::VectorXd *grad, Eigen::MatrixXd *d2)
{
	unsigned num_subpix=__num_x_subpix*__num_y_subpix;
	assert( unsigned( s.size() ) == num_subpix );
	calculate_ws_inverse_f(s);
	Eigen::ArrayXXd xif=__xi.array()*__f.array();
	Eigen::VectorXd fbar=xif.rowwise().sum();
	Eigen::ArrayXXd fbar_minus_f=-__f.array();
	unsigned max_s=__w.size(), max_t=__w[0].size();
	for(unsigned s=0; s<max_s; s++) fbar_minus_f.row(s)+=fbar(s);

    if ( variance ) {
        *variance=__q*(__xi.array() *
            fbar_minus_f.square()).rowwise().sum().matrix();
    }
            
/*		*variance=__q*((__xi.array()*__f.array().square()).rowwise().sum()-
				xif.rowwise().sum().square()).matrix();*/
	if(!grad && !d2) return;
	if(grad) grad->resize(num_subpix);
	if(d2) d2->resize(num_subpix, num_subpix);

	std::valarray<Eigen::ArrayXXd> lambda(num_subpix);
	std::valarray<Eigen::ArrayXd> lambda_bar(num_subpix);
	calculate_lambda(lambda, lambda_bar, num_subpix);
	for(unsigned k=0; k<num_subpix; k++) {
		if(grad)
			(*grad)(k)=2*__q*(__xi.array()*lambda[k]*fbar_minus_f
					).rowwise().sum().matrix();
		if(d2==NULL) continue;
		for(unsigned l=0; l<=k; l++) {
			Eigen::ArrayXXd omega(max_s, max_t);
            calculate_omega( omega, k, l );
			(*d2)(k,l)=-2.0*__q*((__xi.array()*(
					omega*fbar_minus_f-lambda[k]*lambda[l])).rowwise().sum()
					+ lambda_bar[k]*lambda_bar[l]).matrix();
			if(l<k) (*d2)(l,k)=(*d2)(k,l);
		}
	}
}

#endif
