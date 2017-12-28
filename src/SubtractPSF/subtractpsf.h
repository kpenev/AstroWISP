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

class LIB_PUBLIC SubtractPSFCfg : public CommandLineConfig {
private:
	///Describes the available command line options.
	void describe_options();

	///The part of the help describing the usage and purpose (no options).
	std::string usage_help(const std::string &prog_name) const;
public:
	SubtractPSFCfg(
			///The number of arguments on the command line
			///(+1 for the executable)
			int argc,
			
			///A C style array of the actual command line arguments.
			char **argv)
	{parse(argc, argv);}
};

#endif /*SUBTRACT_PSF_H__*/
