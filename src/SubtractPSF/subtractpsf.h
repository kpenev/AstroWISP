#ifndef SUBTRACT_PSF_H__
#define SUBTRACT_PSF_H__

#include "../Core/SharedLibraryExportMacros.h"
#include "CommandLineUtil.h"
#include "H5IODataTree.h"
#include "../IO/SubPixHDF5File.h"
#include "../PSF/DataTreeCalculations.h"
#include "FitsImage.h"
#include "H5Cpp.h"
#include "PiecewiseBicubicPSFMap.h"

///Command line parser for subtracting PSF models from images.
class LIB_PUBLIC SubtractPSFCfg : public CommandLineConfig {
private:
	///Describes the available command line options.
	void describe_options();

	///The part of the help describing the usage and purpose (no options).
	std::string usage_help(
        ///The name with which this tool was invoked.
        const std::string &prog_name
    ) const;
public:
    ///Parse the command line to this object.
	SubtractPSFCfg(
        ///The number of arguments on the command line
        ///(+1 for the executable)
        int argc,

        ///A C style array of the actual command line arguments.
        char **argv
    )
	{parse(argc, argv);}
};

#endif /*SUBTRACT_PSF_H__*/
