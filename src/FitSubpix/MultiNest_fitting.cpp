/**\file
 *
 * \brief The definitions of some of the functions and methods used by the
 * MultiNest based fitting for the sub-pixel map.
 *
 * \ingroup FitSubpix
 */

#include "MultiNest_fitting.h"
#include "MultiNest_v2.14/example_eggbox_C++/multinest.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

namespace MultiNestFit {
	///The subpixel map.
	MultiNestSubPixelMap subpix_map;

	///Calculates the sum of variances for all sources.
	Variance var;

	void log_likelihood(
			double* Cube,
			int&    ndim,
			int&    /*npars*/,
			double& lnew
			)
	{
		Eigen::VectorXd sensitivities=x_to_s(
				Eigen::Map<const Eigen::VectorXd>(Cube, ndim));
		for( int i = 0; i < sensitivities.size(); i++ ) {
			Cube[i] = sensitivities( i );
		}

		double var_value=var.variance(sensitivities);
		std::cout.precision(16);
		std::cout << std::setw(25) << var_value << std::endl;
		lnew=-std::log(var_value+1.0);
	}

	void MultiNestSubPixelMap::set_sensitivities(const double *MNCube)
	{
		unsigned long x_res=x_resolution(), y_res=y_resolution(), 
		nsubpix=x_res*y_res;
		Eigen::VectorXd s=x_to_s(Eigen::Map<const Eigen::VectorXd>(MNCube, 
					nsubpix-1));
		for(unsigned i=0; i<nsubpix; i++) (*this)(i%x_res, i/x_res)=s(i);
	}

	///Outputs the given sensitivities to stdout.
	void dump_sensitivities(double *sensitivities)
	{
		std::ios_base::fmtflags orig_flags=std::cout.flags();
		std::streamsize orig_precision=std::cout.precision();
		std::cout.precision(5);
		std::cout.setf(std::ios::scientific);
		for(unsigned long y=0; y<subpix_map.y_resolution(); y++) {
			for(unsigned long x=0; x<subpix_map.x_resolution()-1; x++)
				std::cout << sensitivities[x+subpix_map.x_resolution()*y]
					<< ' ';
			std::cout << sensitivities[subpix_map.x_resolution()*(y+1)-1]
				<< std::endl;
		}
		std::cout.flags(orig_flags);
		std::cout.precision(orig_precision);
	}

	void dumper(
			int&        nSamples,
			int&        nlive,
			int&        nPar,
			double**    physLive,
			double**    posterior,
			double**    paramConstr,
			double&     /*maxLogLike*/,
			double&     /*logZ*/,
			double&     /*logZerr*/
			)
	{
		std::cout << "Mean sensitivities:" << std::endl;
		dump_sensitivities(paramConstr[0]);
		std::cout << "Standard deviation of sensitivities:" << std::endl;
		dump_sensitivities(paramConstr[0]+nPar);
		std::cout << "Maximum likelihood sensitivities:" << std::endl;
		dump_sensitivities(paramConstr[0]+2*nPar);
		std::cout << "Maximum-a-posteriory sensitivities:" << std::endl;
		dump_sensitivities(paramConstr[0]+3*nPar);
		static int dump_index=0;
		std::ostringstream fname;
		fname << "output/MultiNest_dump_" << std::setfill('0') << std::setw(5)
			<< dump_index << '.';

		std::cout << "Writing active points and likelihood values:" << std::endl;
		std::ofstream active_os((fname.str()+".active").c_str());
		for(int pt_ind=0; pt_ind<nlive; pt_ind++) {
			active_os << "\t";
			for(int par_ind=0; par_ind<nPar; par_ind++)
				active_os << std::setw(25) << physLive[0][par_ind+pt_ind*nPar];
			active_os << std::setw(25) << physLive[0][nlive*nPar+pt_ind]
				<< std::endl;
		}
		active_os.close();
		std::cout << "Writing posterior distribution:" << std::endl;
		std::ofstream post_os((fname.str()+".posterior").c_str());
		for(int pt_ind=0; pt_ind<nSamples; pt_ind++) {
			post_os << "\t";
			for(int par_ind=0; par_ind<nPar; par_ind++)
				post_os << std::setw(25) << posterior[0][par_ind+pt_ind*nPar];
			post_os << std::setw(25) << posterior[0][nSamples*nPar+pt_ind] 
				<< std::setw(25) << posterior[0][nSamples*(nPar+1)+pt_ind]
				<< std::endl;
		}
		post_os.close();
		dump_index++;
	}

	void fit(
			const int&          x_split,
			const int&          y_split,
			const StringList&   frame_filenames,
			const StringList&   source_filenames,
			double              aperture,
			)
	{
		// set the MultiNest sampling parameters
		int mmodal = 1;					//do mode separation?

		int ceff = 1;					//run in constant efficiency mode?

		int nlive = 5000;				//number of live points

		double efr = 0.5;				//set the required efficiency

		double tol = 0.01;				//tol, defines the stopping criteria

		int ndims = x_split*y_split-1;	//dimensionality (no. of free 
		//parameters)

		int nPar = ndims+1;				//total no. of parameters including
		//free & derived parameters

		int nClsPar = ndims;			//no. of parameters to do mode 
		//separation on

		int updInt = 1;					//after how many iterations feedback
		//is required & the output files 
		//should be updated
		//note: posterior files are updated &
		//dumper routine is called after 
		//every updInt*10 iterations

		double Ztol = -1E90;			//all the modes with logZ < Ztol are
		//ignored

		int maxModes = 100;				//expected max no. of modes (used 
		//only for memory allocation)

		int pWrap[ndims];				//which parameters to have periodic
		//boundary conditions?
		for(int i = 0; i < ndims; i++) pWrap[i] = 0;

		char root[100] = "chains/1-";	//root for output files

		int seed = -1;					//random no. generator seed, if < 0
		//then take the seed from system 
		//clock

		int fb = 1;						//need feedback on standard output?

		int resume = 1;					//resume from a previous job?

		int outfile = 0;				//write output files?

		int initMPI = 1;				//initialize MPI routines?, relevant
		//only if compiling with MPI
		//set it to F if you want your main 
		//program to handle MPI 
		//initialization

		double logZero = -1E90;			//points with loglike < logZero will
		//be ignored by MultiNest

		//	int context = 0;				//not required by MultiNest, any
		//additional information user wants 
		//to pass.

		//Initialize the chi squared variable and supixel map.
		subpix_map.set_resolution(x_split, y_split);
		var.initialize(
				frame_filenames,
				source_filenames,
				aperture,
				x_split,
				y_split
				);

		// calling MultiNest
		nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, 
				maxModes, updInt, Ztol, root, seed, pWrap, fb, resume,
				outfile, initMPI, logZero, log_likelihood, 
				dumper);
	}
};
