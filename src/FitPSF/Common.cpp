#include "Common.h"
#include "LinearSource.h"

namespace FitPSF {

    std::ostream &operator<<(std::ostream &os,
                             const SourceDropReason &reason)
    {
        switch(reason) {
            case FEW_PIXELS: os << "FEW_PIXELS"; break;
            case MANY_PIXELS: os << "MANY_PIXELS"; break;
            case TOO_BIG: os << "TOO_BIG"; break;
            case OVERLAP: os << "OVERLAP"; break;
            case NON_POINT_SOURCE: os << "NON_POINT_SOURCE"; break;
            case BAD_BACKGROUND: os << "BAD_BACKGROUND"; break;
            case PAST_MAX_SOURCES: os << "PAST_MAX_SOURCES"; break;
            case MANY_SATURATED: os << "MANY_SATURATED"; break;
            default : assert(false);
        }
        return os;
    }

    void add_new_source(
            Image<LinearSource>                     &image,
            const Core::SubPixelMap                 *subpix_map,
            const PSF::PiecewiseBicubic             &psf,
            double                                   alpha,
            double                                   max_circular_aperture,
            const std::string                       &output_fname,
            bool                                     cover_psf,
            const Core::SourceLocation              &location,
            const Eigen::VectorXd                   &psf_terms,
            const Background::Source                &srcbg,
            size_t                                   source_assignment_id,
            LinearSourceList                        &destination
    )
    {
        if(cover_psf) {
            // assing pixels overlapping with the PSF grid
            destination.push_back(
                new LinearSource(psf,
                                 location.id(),
                                 location.x(),
                                 location.y(),
                                 srcbg,
                                 image,
                                 source_assignment_id,
                                 subpix_map,
                                 psf.min_x(),
                                 psf.max_x(),
                                 psf.min_y(),
                                 psf.max_y(),
                                 output_fname)
            );
        } else {
            // assign pixels based on signal to noise
            destination.push_back(
                new LinearSource(psf,
                                 location.id(),
                                 location.x(),
                                 location.y(),
                                 srcbg,
                                 image,
                                 alpha,
                                 source_assignment_id,
                                 subpix_map,
                                 max_circular_aperture,
                                 output_fname)
            );
        }
        Eigen::VectorXd &source_expansion_terms =
            destination.back()->expansion_terms();
        source_expansion_terms.resize(psf_terms.size());
        source_expansion_terms = psf_terms;
    }
} //End FitPSF namespace.
