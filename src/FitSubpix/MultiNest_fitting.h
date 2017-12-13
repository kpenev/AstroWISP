/**\file
 *
 * \brief Declarations of the classes and functions related to the MultiNest
 * based fitting for the sub-pixel map.
 *
 * \ingroup FitSubpix
 */

#ifndef __MULTI_NEST_FIT_H
#define __MULTI_NEST_FIT_H

#include "SubPixelMap.h"
#include "FittingUtil.h"

/**\brief The namespace isolating anything exclusively needed by the
 * MultiNest fitting for the sub-pixel map.
 *
 * \ingroup FitSubpix
 */
namespace MultiNestFit {
	/**\brief A subpixel map that sets its sensitivities from a MultiNest
	 * unit cube.
	 *
	 * \ingroup FitSubpix
	 */
	class MultiNestSubPixelMap : public SubPixelMap {
	public:
		///Just name the map, resultion must be set later.
		MultiNestSubPixelMap() : SubPixelMap("MultiNest sub pixel map") {}

		///Construct a map with the given resolution.
		MultiNestSubPixelMap(unsigned long x_split, unsigned long y_split) :
			SubPixelMap(x_split, y_split, "MultiNest sub pixel map") {}

		///\brief Sets the sensitivities from the given MultiNest unit cube.
		///
		///The dimensionality of the cube should be one less than the number
		///of subpixels. The mapping of the unit cube to subpixel
		///sensitivities is such that the prior probability for all subpixel
		///maps is the same. If the values of the unit cube are u_i then the
		///sensitivities s_i are calculated as follows: 
		///s_i=c_i{1-(1-u_i)^[1/(N-i)]}, 
		///with c_i=1-SUM(s_j, j<i), except for s_N which is set to c_N. With
		///the sensitivities ordered such that x changes faster with i.
		void set_sensitivities(const double *MNCube);
	};


	///The actual function to be maximized by MultiNest: -Log(Chi2).
	void log_likelihood(
				///\brief On input: the ndim parameters in unit-hypercube
				///On output: all sensitivities, including the last one.
				double *Cube,

				///One less than the number of supixels, input only.
				int &ndim, 

				///\brief Total number of subpixels (number of output
				///parameters), input only.
				int &npars, 

				///\brief The value of the function to maximize for the given
				///hypercube values.
				double &lnew);

	///Perform the actual fitting.
    void fit(
        const int&          x_split,
        const int&          y_split,
        const StringList&   frame_filenames,
        const StringList&   source_filenames,
        double              aperture,
    );

	///\brief The dumper routine will be called every updInt*10 iterations.
	///
	///MultiNest doesn not need to the user to do anything. User can use the
	///arguments in whichever way he/she wants
	void dumper(
		///Total number of samples in posterior distribution.
		int &nSamples, 

		///total number of live points
		int &nlive, 

		///Total number of parameters (free + derived).
		int &nPar, 
		
		///2D array containing the last set of live points (physical
		///parameters plus derived parameters) along with their loglikelihood
		///values: physLive[1][nlive * (nPar + 1)]
		double **physLive, 

		///posterior distribution containing nSamples points. Each sample has
		///nPar parameters (physical + derived) along with the their loglike
		///value & posterior probability: posterior[1][nSamples * (nPar + 2)]
		double **posterior, 

		///The contents is:
		/// - paramConstr[0][0] to paramConstr[0][nPar - 1]: mean values of 
		///   the parameters
		/// - paramConstr[0][nPar] to paramConstr[0][2*nPar - 1]: standard
		///   deviation of the parameters
		/// - paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1]: best-fit
		///   (maxlike) parameters
		/// - paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1]: MAP
		///   (maximum-a-posteriori) parameters
		double **paramConstr, 

		///maximum loglikelihood value
		double &maxLogLike, 
		
		///log evidence value
		double &logZ, 

		///error on log evidence value
		double &logZerr);
}

#endif
