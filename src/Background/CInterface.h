/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * Background library.
 *
 * \ingroup Background
 */

#include "MeasureAnnulus.h"
#include "../Core/Image.h"
#include "../Core/CInterface.h"

extern "C" {
    ///Opaque struct to cast to/from Background::MeasureAnnulus.
    struct BackgroundMeasureAnnulus;

    ///\brief Create an object for measuring the background in an annulus around
    ///each source.
    BackgroundMeasureAnnulus *create_background_extractor(
        ///Inner radius of the annulus used for background determination.
        double inner_radius,

        ///Outer radius of the annulus used for background determination.
        double outer_radius, 

        ///Size of the area to exclude from other source's annuli.
        double exclude_aperture,

        ///The image to determine the backgrounds on.
        const Image *image, 

        ///The confidence to require of the error estimate.
        double error_confidence
    );

    ///Destroy a previously created background extractor.
    void destroy_background_extractor(
        ///The background extractor to destroy. Must have been created using
        ///create_background_extractor().
        BackgroundMeasureAnnulus *extractor
    );

    ///\brief Register another source around which to exclude pixels from
    ///background estimation.
    void add_source_to_background_extractor(
        ///The background extractor to add the source to. Must have been created
        ///using create_background_extractor().
        BackgroundMeasureAnnulus *extractor,

        ///The x coordinate (within the image) of the source to add.
        double x,

        ///The y coordinate (within the image) of the source to add.
        double y
    );

    ///Estimate the background at a specified position.
    void measure_background(
        ///The background extractor to use.
        const BackgroundMeasureAnnulus *extractor,

        ///The x coordinate (within the image) for which to estimate the
        ///background.
        double x,

        ///The y coordinate (within the image) for which to estimate the
        ///background.
        double y,

        ///Overwritten with the best estimate of the background at (x, y).
        double *value,

        ///Overwritten with the best estimate of the uncertainty in the
        ///backgronud determination.
        double *error,

        ///Overwritten with the number of pixels that contributed to the
        ///determination of the value and error.
        unsigned *pixels
    );

    ///\brief (Re)start iterating over the sources currently registered with an
    ///extractor measuring each one's background.
    void restart_background_iteration(
        ///The background extractor to use.
        BackgroundMeasureAnnulus *extractor
    );

    ///\brief Get the background under the current source and advance to the
    ///next source.
    ///
    ///Return valuen indicates whether there are more sources pending.
    bool get_next_background(
        ///The background extractor to use.
        BackgroundMeasureAnnulus *extractor,

        ///Overwritten with the best estimate of the background at (x, y).
        double *value,

        ///Overwritten with the best estimate of the uncertainty in the
        ///backgronud determination.
        double *error,

        ///Overwritten with the number of pixels that contributed to the
        ///determination of the value and error.
        unsigned *pixels
    );
};//End extern "C".
