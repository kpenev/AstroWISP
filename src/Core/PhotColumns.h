/**\file
 *
 * \brief Defines constants to denote the various columns that can be input
 * or output by the various tools and some utility functions.
 *
 * \ingroup Core
 */

#ifndef __PHOT_COLUMNS_H
#define __PHOT_COLUMNS_H

#include "../Core/SharedLibraryExportMacros.h"
#include <cmath>
#include <vector>
#include <string>

///\brief A namespace to isolate the column constants since they are given
///quite common names.
///
///\ingroup SubPixPhot
///\ingroup FitSubpix
///\ingroup FitPSF
namespace Phot {
	///\brief constants to denote the various columns that can be input or
	///output by the SubPixPhot tool
	enum LIB_PUBLIC Columns {
		id,			///< The HAT-id.
		x,			///< The x coordinate in the input fits image.
		y,			///< The y coordinate in the input fits image.
		S,			///< The S paramater of an elliptical gaussian PSF.
		D,			///< The D paramater of an elliptical gaussian PSF.
		K,			///< The K paramater of an elliptical gaussian PSF.
		A,			///< The amplitude of an elliptical gaussian PSF.
		bg,			///< The background under the source.
		bg_err,		///< The uncertainty of the background under the source.
		flux,		///< The flux of a source.
		flux_err,	///< The uncertainty in the flux of a source.
		mag,		///< The magnitude of a source.
		mag_err,	///< The uncertainty in the magnitude of a source.
		chi2,		///< \f$\chi^2\f$ of a PSF fit to source.
		sn,			///< Signal to noise of the source.
		FIRST_INT_COLUMN,
		npix=FIRST_INT_COLUMN,///< The number of pixels assigned to a source.
        nbgpix,    ///< The number of background pixels assigned to a source.
		flag,		///< Flag indicating the quality of the flux measurement.
		enabled,	///< A switch to enable/disable the source.

#ifdef DEBUG
		time,
#endif
		unknown,	///< A column that is not used but is in a file.

		///How many different types of columns are there.
		num_recognized_columns
	};

    ///\brief Names to use for printing Columns to streams.
	class LIB_LOCAL ColumnNamesVector : public std::vector<std::string> {
	public:
		ColumnNamesVector() :
			std::vector<std::string>(num_recognized_columns)
		{
			std::vector<std::string> &vector=*this;
			vector[id]="ID";
			vector[x]="x";
			vector[y]="y";
			vector[S]="S";
			vector[D]="D";
			vector[K]="K";
			vector[A]="Amplitude";
			vector[bg]="Background";
			vector[bg_err]="BackgroundErr";
			vector[flux]="Flux";
			vector[flux_err]="FluxErr";
			vector[mag]="Magnitude";
			vector[mag_err]="MagnitudeErr";
			vector[chi2]="Chi2";
			vector[sn]="SignalToNoise";
			vector[npix]="NPix";
			vector[nbgpix]="BackgroundNPix";
			vector[flag]="QualityFlag";
			vector[enabled]="Enabled";
			vector[unknown]="Unknown";
#ifdef DEBUG
			vector[time]="ProcessingTime";
#endif
		}
	};

    ///Instantiate of ColumnNamesVector so it can be used.
	const ColumnNamesVector column_name;
}

///Calculate the magnitude corresponding to a given flux.
inline LIB_LOCAL double magnitude(
    ///The flux to convert to a magnitude.
    double flux,

    ///The magnitude corresponding to a flux of 1.
    double mag_1adu
)
{
	return mag_1adu-2.5*log10(flux);
}

///The error in the magnitude corresponding to a given flux.
inline LIB_LOCAL double magnitude_error(
    ///The flux being converted to a magnitude.
    double flux,

    ///The error estimate of the flux.
    double flux_error
)
{
	return -2.5*log10(1.0-flux_error/flux);
}

#endif
