/**\file
 *
 * \brief Declares the command line parser for the VisualizePSFMap tool.
 */

#ifndef __DISPLAY_PSF_MAP_H
#define __DISPLAY_PSF_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "FitsImage.h"
#include "PiecewiseBicubicPSFMap.h"
#include "EllipticalGaussianPSFMap.h
#include "H5IODataTree.h"
#include "CommandLineUtil.h"
#include "Typedefs.h"
#include "Error.h"
#include <list>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdio>
#include <fstream>

///Command line parser for Visualizing PSF maps.
class LIB_PUBLIC VisualizePSFMapConfig : public CommandLineConfig {
private:
	///Describes the available command line options.
	void describe_options();

public:
	///Parse the command line.
	VisualizePSFMapConfig(
			///The number of arguments on the command line
			///(+1 for the executable)
			int argc,

			///A C style array of the actual command line arguments.
			char **argv)
	{
        parse(argc, argv);
        PSF::EllipticalGaussian::set_default_precision(
            operator[]("psf.sdk.rel-int-precision").as<double>(),
            operator[]("psf.sdk.abs-int-precision").as<double>()
        );
        PSF::EllipticalGaussian::set_default_max_exp_coef(
            operator[]("psf.sdk.max-exp-coef").as<double>()
        );
    }
};

#endif
