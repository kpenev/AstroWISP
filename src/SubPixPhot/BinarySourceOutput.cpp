/**\file
 *
 * \brief Define some of the methods of the BinarySourceOutput class.
 *
 * \ingroup SubPixPhot
 */

#include "BinarySourceOutput.h"

namespace SubPixPhot {

    void BinarySourceOutput::set_required_columns()
    {
        __required = false;
        unsigned num_flux = 0,
                 num_flux_err = 0,
                 num_mag = 0,
                 num_mag_err = 0,
                 num_flag = 0;
        for(
            std::list<Phot::Columns>::const_iterator ci = __columns->begin();
            ci != __columns->end();
            ++ci
        ) {
            __required[*ci] = true;
            if(*ci == Phot::flux) num_flux++;
            else if(*ci == Phot::flux_err) ++num_flux_err;
            else if(*ci == Phot::mag) ++num_mag;
            else if(*ci == Phot::mag_err) ++num_mag_err;
            else if(*ci == Phot::flag) ++num_flag;
        }
        if(__required[Phot::mag]) __required[Phot::flux] = true;
        if(__required[Phot::mag_err]) 
            __required[Phot::flux] = __required[Phot::flux_err] = true;
        __num_apertures = std::max(
            num_flux, 
            std::max(
                num_flux_err, 
                std::max(
                    num_mag,
                    std::max(num_mag_err, num_flag)
                )
            )
        );
    }

    std::string BinarySourceOutput::aperture_indep_format()
    {
        std::string result;
        for(
            std::list<Phot::Columns>::const_iterator ci = __columns->begin();
            ci != __columns->end();
            ++ci
        ) {
            switch(*ci) {
                case Phot::id : result += 'I'; break;
                case Phot::x : result += 'X'; break;
                case Phot::y : result += 'Y'; break;
                case Phot::S :
                case Phot::D :
                case Phot::K :
                case Phot::A :
                case Phot::enabled :
                case Phot::chi2 :
                case Phot::npix :
                case Phot::unknown :
                           throw Error::IO("Not all requested output "
                                           "columns are supported.");
                default : ;
            }
        }
        return result;
    }

    std::string BinarySourceOutput::aperture_dep_format()
    {
        std::string result;
        for(
            std::list<Phot::Columns>::const_iterator ci = __columns->begin();
            ci != __columns->end();
            ci++
        ) {
            switch(*ci) {
                case Phot::flux : result += 'F'; break;
                case Phot::flux_err : result += 'f'; break;
                case Phot::mag : result += 'M'; break;
                case Phot::mag_err : result += 'm'; break;
                case Phot::flag : result += 'S'; break;
                default :;
            }
        }
        return result;
    }

    BinarySourceOutput::BinarySourceOutput(
        IO::binostream &outstream,
        double mag_1ADU,
        const std::list<Phot::Columns> &output_columns
    ) :
        __outstream(&outstream),
        __required(Phot::num_recognized_columns),
        __precision(Phot::num_recognized_columns),
        __mag_1ADU(mag_1ADU)
    {
        set_columns(output_columns);
        __precision[Phot::x] =
            __precision[Phot::y] =
            __precision[Phot::S] =
            __precision[Phot::D] =
            __precision[Phot::K] = 2;
        __precision[Phot::A] = 0;
        __precision[Phot::bg] = __precision[Phot::bg_err] = 2;
        __precision[Phot::flux] = __precision[Phot::flux_err] = 0;
        __precision[Phot::mag] = __precision[Phot::mag_err] = 5;
    }

    void BinarySourceOutput::set_columns(
        const std::list<Phot::Columns> &output_columns
    )
    {
        __columns = &output_columns;
        set_required_columns();
    }

    void BinarySourceOutput::operator()(
        const std::list<IO::OutputSDKSource> &sources
    )
    {
        size_t num_sources = sources.size();
        std::valarray<unsigned long>
                id_field(__required[Phot::id] ? num_sources : 0),
                id_source(__required[Phot::id] ? num_sources : 0);
        std::valarray<double>
            x(__required[Phot::x] ? num_sources : 0), 
            y(__required[Phot::y] ? num_sources : 0),
            s(__required[Phot::S] ? num_sources : 0),
            d(__required[Phot::D] ? num_sources : 0),
            k(__required[Phot::K] ? num_sources : 0),
            amp(__required[Phot::A] ? num_sources : 0),
            bg(__required[Phot::bg] ? num_sources : 0),
            bg_err(__required[Phot::bg_err] ? num_sources : 0);

        std::valarray< std::valarray<double> > 
            flux(
                std::valarray<double>(
                    __required[Phot::flux] ? num_sources : 0
                ),
                __num_apertures
            ),
            flux_err(
                std::valarray<double>(
                    __required[Phot::flux_err] ?  num_sources : 0
                ),
                __num_apertures
            ),
            mag(__num_apertures),
            mag_err(__num_apertures);
        std::valarray< std::valarray<char> > flag(
            std::valarray<char>(__required[Phot::flag] ? num_sources : 0),
            __num_apertures
        );
        size_t index = 0;
        for(
            std::list<IO::OutputSDKSource>::const_iterator
                si = sources.begin();
            si != sources.end();
            si++
        ) {
            if(__required[Phot::id]) {
                if ( si->id().is_hatid() ) {
                    id_field[index] = si->id().field();
                    id_source[index] = si->id().source();
                }
                else {
                    // TODO: invalid HAT ID: add warning?
                    id_field[index] = 0;
                    id_source[index] = index;
                }
            }
            if(__required[Phot::x])
                x[index] = si->x();
            if(__required[Phot::y])
                y[index] = si->y();
            if(__required[Phot::S])
                s[index] = si->psf_s();
            if(__required[Phot::D])
                d[index] = si->psf_d();
            if(__required[Phot::K])
                k[index] = si->psf_k();
            if(__required[Phot::A])
                amp[index] = si->psf_amplitude();
            if(__required[Phot::bg])
                bg[index] = si->background().value();
            if(__required[Phot::bg_err])
                bg_err[index] = si->background().error();
            for(unsigned ap_ind = 0; ap_ind < __num_apertures; ++ap_ind) {
                if(__required[Phot::flux])
                    flux[ap_ind][index] = si->flux()[ap_ind].value();
                if(__required[Phot::flux_err])
                    flux_err[ap_ind][index] = si->flux()[ap_ind].error();
                if(__required[Phot::flag])
                    flag[ap_ind][index] = si->flux()[ap_ind].flag();
            }
            ++index;
        }
        for(unsigned ap_ind = 0; ap_ind < __num_apertures; ++ap_ind) {
            if(__required[Phot::mag]) {
                mag[ap_ind].resize(flux.size());
                mag[ap_ind] = __mag_1ADU - 2.5 * log10(flux[ap_ind]);
            }
            if(__required[Phot::mag_err]) {
                mag_err.resize(flux.size());
                mag_err[ap_ind] = -2.5 * log10(
                    1.0 - flux_err[ap_ind] / flux[ap_ind]
                );
            }
        }
        unsigned flux_ind = 0,
                 flux_err_ind = 0,
                 flag_ind = 0,
                 mag_ind = 0,
                 mag_err_ind = 0;
        *__outstream << aperture_indep_format().c_str() << '\0'
                     << aperture_dep_format().c_str() << '\0'
                     << char(1)
                     << static_cast<char>(__num_apertures);
        for(
            std::list<Phot::Columns>::const_iterator ci = __columns->begin();
            ci != __columns->end();
            ++ci
        ) {
            __outstream->precision(__precision[*ci]);
            switch(*ci) {
                case Phot::id : *__outstream << id_field << id_source; break;
                case Phot::x : *__outstream << x; break;
                case Phot::y : *__outstream << y; break;
                case Phot::S : *__outstream << s; break;
                case Phot::D : *__outstream << d; break;
                case Phot::K : *__outstream << k; break;
                case Phot::A : *__outstream << amp; break;
                case Phot::bg : *__outstream << bg; break;
                case Phot::bg_err : *__outstream << bg_err; break;
                case Phot::flux : *__outstream << flux[flux_ind++]; break;
                case Phot::flux_err : 
                                  *__outstream << flux_err[flux_err_ind++];
                                  break;
                case Phot::mag : *__outstream << mag[mag_ind++]; break;
                case Phot::mag_err :
                                 *__outstream << mag_err[mag_err_ind++];
                                 break;
                case Phot::flag : *__outstream << flag[flag_ind++]; break;
                default :;

            };
        }
    }

} //End SubPixPhot namespace.
