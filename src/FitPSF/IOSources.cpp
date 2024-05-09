#include "IOSources.h"

namespace FitPSF {
    IOSources::IOSources(const char *fits_fname,
                         char **source_ids,
                         const double *source_coordinates,
                         const double *psf_terms,
                         const bool *enabled,
                         unsigned long num_sources,
                         unsigned long num_terms) :
        __locations(num_sources),
        __psf_terms(num_terms, num_sources),
        __fits_fname(fits_fname),
        __output_fname(fits_fname)
    {
        __psf_terms = Eigen::Map<Eigen::MatrixXd>(
            const_cast<double*>(psf_terms),
            num_terms,
            num_sources
        );

        if(enabled)
            __enabled.assign(enabled, enabled + num_sources);

        for(
            unsigned long source_ind = 0;
            source_ind < num_sources;
            ++source_ind
        ) {
            __locations[source_ind].id() = Core::SourceID(
                source_ids[source_ind]
            );
             __locations[source_ind].x() = source_coordinates[2 * source_ind];
             __locations[source_ind].y() = source_coordinates[2 * source_ind
                                                              +
                                                              1];
        }
    }
}
