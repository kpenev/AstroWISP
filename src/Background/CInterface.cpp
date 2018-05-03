/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup Background
 */

#include "CInterface.h"

BackgroundMeasureAnnulus *create_background_extractor(
    double inner_radius,
    double outer_radius, 
    double exclude_aperture,
    const CoreImage *image,
    double error_confidence
)
{
    return reinterpret_cast<BackgroundMeasureAnnulus*>(
        new Background::MeasureAnnulus(
            inner_radius,
            outer_radius,
            exclude_aperture,
            *reinterpret_cast<const Core::Image<double> *>(image),
            error_confidence
        )
    );
}

void destroy_background_extractor(BackgroundMeasureAnnulus *extractor)
{
    delete reinterpret_cast<Background::MeasureAnnulus*>(extractor);
}

void add_source_to_background_extractor(BackgroundMeasureAnnulus *extractor,
                                        double x,
                                        double y)
{
    reinterpret_cast<Background::MeasureAnnulus*>(extractor)->add_source(x, y);
}

void measure_background(BackgroundMeasureAnnulus *extractor,
                        double x,
                        double y,
                        double *value,
                        double *error,
                        unsigned *pixels)
{
    Background::Source source =
        reinterpret_cast<Background::MeasureAnnulus*>(extractor)->operator()(x,
                                                                             y);
    *value = source.value();
    *error = source.error();
    *pixels = source.pixels();
}

void restart_background_iteration(BackgroundMeasureAnnulus *extractor)
{
    reinterpret_cast<Background::MeasureAnnulus*>(
        extractor
    )->jump_to_first_source();
}

bool get_next_background(BackgroundMeasureAnnulus *extractor,
                         double *value,
                         double *error,
                         unsigned *pixels)
{
    Background::MeasureAnnulus *real_extractor =
        reinterpret_cast<Background::MeasureAnnulus*>(extractor);
    Background::Source source = real_extractor->operator()();

    *value = source.value();
    *error = source.error();
    *pixels = source.pixels();
    return real_extractor->next_source();
}
