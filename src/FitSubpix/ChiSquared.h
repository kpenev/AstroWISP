#ifndef __CHI_SQUARED_H
#define __CHI_SQUARED_H

/**\addtogroup FitSubpix
 * @{
 */

/**\file
 * \brief Contains the declaration of the ChiSquared class.
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "FitsImage.h"
#include "Source.h"
#include "PhotColumns.h"
#include "CommandLineUtil.h"
#include "SubPixPhotIO.h"
#include "SubPixelMap.h"
#include <vector>
#include <list>
#include <valarray>
#include <string>

///\brief The function minimized when fitting for the sub-pixel sensitivities
///(simplex method only).
///
///This class is almost obsolete. It is only used by the simplex fitting
///method, which itself is obsolete.
class LIB_LOCAL ChiSquared {
private:
	///The fits frames to use for photometry.
	std::vector< FitsImage<double> > __images;

    ///The information read from the input files.
    std::vector<H5IODataTree> __data_trees;

	mutable std::valarray< std::valarray<double> > 
        ///Set of flux measurements for each source common to all images.
        __flux_values, 

        ///Error in flux measurements for each source common to all images.
        __flux_errors; 
	
	///The aperture to use when doing photometry and when masking noughboring
	///sources.
	double __aperture;

	///Check that a data tree contains the same sources as the first one.
	void verify_sources(
            ///The data tree whose source list should be compared to
            ///__data_trees[0]
            const H5IODataTree &this_tree,

            ///The name of the file from which this tree was read (only used
            ///in the message generated if a mismatch in the source list is
            ///found).
            const std::string &psf_fname);

    ///\brief Read data from a psf file to __data_trees.
    void read_psf_file(
            ///The name of the file.
            const std::string &filename,
            
            ///The index witihn __data_trees to fill with the new data.
            unsigned image_index);

    ///\brief Measure fluxes of all image sources and add to correspnoding
    ///data tree.
    void measure_image_fluxes(
            ///The index of the image to process.
            unsigned image_index,

            ///The sub-pixel map to use (assumed properly normalized).
            const SubPixelMap &subpix_map,

            ///A list containing only one number - the aperture for which
            ///ChiSquared is to be calculated.
            const std::list<double> &aperture_list);
protected:
	///Calculates the variance (appropriately scaling all deviations by the
	///error).
	virtual double chi2() const;

public:
	///Postpone proper initialization until later.
	ChiSquared() {};

	///Create a fully functional Chi2.
	ChiSquared(
			///Frames to use for fitting the subpixel structure.
			const StringList &frame_filenames,

			///Source lists: one for each frame.
			const StringList &psf_filenames,

			///The aperture to use when deriving the flux.
			double aperture)
	{
        initialize(frame_filenames, 
                   psf_filenames,
                   aperture);
    }


	///Delayed initialization or reset. Same arguments as the constructor.
	void initialize(const StringList &frame_filenames,
					const StringList &psf_filenames,
					double aperture);

	///The value of the normalized variance for the given subpixel
	///sensitivity map
	double operator()(const SubPixelMap &subpix_map) const;
};

/** @} */

#endif
