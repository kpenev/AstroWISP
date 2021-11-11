/**\file
 *
 * \brief Declare a class for extracting source background from an annulus
 * around each source excluding areas around nearby sources.
 *
 * \ingroup Background
 */

#ifndef __MEASURE_ANNULUS_H
#define __MEASURE_ANNULUS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "MeasureExcludingSources.h"
#include <valarray>
#include <cmath>

namespace Background {

    ///\ingroup SubPixPhot
    ///\ingroup FitSubpix
    ///\brief Measure the background in an annulus around each source.
    ///
    ///Uses the pixels whose centers lie between two circles centered on each
    ///source to estimate the background (median) and its error (quantiles in
    ///the individual pixel values).
    ///
    ///Pixels within the annulus for a given source are not included in the
    ///background determination if they fall inside the exclusion radius of
    ///another source.
    class LIB_PUBLIC MeasureAnnulus : public MeasureExcludingSources {
    private:
        double
            ///The inner radius of the region around a source used to
            ///deterime the background.
            __inner_radius,

            ///The outer radius of the region around a source used to
            ///deterime the background.
            __outer_radius;

        ///What fraction of the pixels on which the background is based should
        ///fall within the error range.
        double __error_confidence;

    public:
        ///\brief Use the given radii to determive background.
        MeasureAnnulus(
            ///Inner radius of the annulus used for background determination.
            double inner_radius,

            ///Outer radius of the annulus used for background determination.
            double outer_radius,

            ///Size of the area to exclude from other source's annuli.
            double exclude_aperture,

            ///The image to determine the backgrounds on.
            const Core::Image<double> &image,

            ///The confidence to require of the error estimate.
            double error_confidence=0.68
        ) :
            MeasureExcludingSources(image, exclude_aperture),
            __inner_radius(inner_radius), __outer_radius(outer_radius),
            __error_confidence(error_confidence)
        {}

        ///Use the given radii to determive background.
        template<class POINT_TYPE>
        MeasureAnnulus(
            ///Inner radius of the annulus used for background determination.
            double inner_radius,

            ///Outer radius of the annulus used for background determination.
            double outer_radius,

            ///Size of the area to exclude from other source's annuli.
            double exclude_aperture,

            ///The image to determine the backgrounds on.
            const Core::Image<double> &image,

            ///The sources in the image.
            const std::list< POINT_TYPE > &sources,

            ///The confidence to require of the error estimate.
            double error_confidence=0.68
        ) :
            MeasureExcludingSources(image, exclude_aperture, sources),
            __inner_radius(inner_radius),
            __outer_radius(outer_radius),
            __error_confidence(error_confidence)
        {}

        ///\brief measure the background for a source at (x,y).
        ///
        ///Estimate the background around the given position as the median of
        ///the points whose centers lie within the annulus and not within the
        ///inner edge of another source's annulus.
        ///
        ///The uncertainty of the determined value is estimated as
        /// \f$\sqrt{\frac{\pi}{2(N-1)}}\f$ times half the
        ///size of the interval centered on the determined background value
        ///that contains the fraction of the points specified at construction
        ///as the error_confidence argument.
        Source operator()(double x, double y) const;

        using MeasureExcludingSources::operator();
    }; //End MeasureAnnulus class.

} //End Background namespace

#endif
