/**\file
 *
 * \brief Define some of the methods of the SourceInput class.
 *
 * \ingroup SubPixPhot
 */

#include "SourceInput.h"

namespace SubPixPhot {

    void SourceInput::set_columns(const std::list<Phot::Columns> &columns)
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Setting "
                  << columns.size()
                  << " source columns." << std::endl;
#endif
        unsigned col=0;
        __id_col = 
            __x_col =
            __y_col =
            __s_col =
            __d_col =
            __k_col =
            __amp_col =
            __bg_col =
            __on_col = -1;
        for(
            std::list<Phot::Columns>::const_iterator ci = columns.begin(); 
            ci!=columns.end();
            ++ci
        ) {
#ifdef TRACK_PROGRESS
            std::cerr << "Found " << Phot::column_name[*ci] << std::endl;
#endif
            switch(*ci) {
                case Phot::id : __id_col = col; break;
                case Phot::x : __x_col = col; break;
                case Phot::y : __y_col = col; break;
                case Phot::S : __s_col = col; break;
                case Phot::D : __d_col = col; break;
                case Phot::K : __k_col = col; break;
                case Phot::A : __amp_col = col; break;
                case Phot::bg : __bg_col = col; break;
                case Phot::enabled : __on_col = col; break;
                default:;
            }
            ++col;
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Set all columns." << std::endl;
#endif
    }

    ///Reads a source from the input stream to src.
    template<class SOURCE_TYPE>
        bool SourceInput::read_source(SOURCE_TYPE &source,
                                      double &s,
                                      double &d,
                                      double &k,
                                      double &amp,
                                      double &bg)
        {
            std::string line("#") ;
            while(line[0] == '#') getline(*__instream, line);
            std::istringstream line_stream(line);
            if(!(*__instream)) {
                __eof = __instream->eof();
                __good = __eof;
                return false;
            }
            double x, y;
            std::string dummy;
            int on = 1;
            Core::SourceID id;
            for(int col = 0; col < __col_num; ++col) {
                if(!line_stream) {
                    __good = false;
                    __eof = __instream->eof();
                    return false;
                }
                if(col == __id_col) line_stream >> id;
                else if(col == __x_col) line_stream >> x;
                else if(col == __y_col) line_stream >> y;
                else if(col == __s_col) line_stream >> s;
                else if(col == __d_col) line_stream >> d;
                else if(col == __k_col) line_stream >> k;
                else if(col == __amp_col) line_stream >> amp;
                else if(col == __bg_col) line_stream >> bg;
                else if(col == __on_col) line_stream >> on;
                else line_stream >> dummy;
            }
            source.id() = id;
            source.x() = x;
            source.y() = y;
            return on;
        }

    ///Read a source with an Elliptical Gaussian PSF.
    SourceInput &SourceInput::operator>>(Core::SDKSource &source)
    {
        int orig_on_col = __on_col; 
        __on_col = -1; 
        double s, d, k, amp = 1.0, bg = 0.0;
        read_source(source, s, d, k, amp, bg); 
        __on_col = orig_on_col;
        source.set_psf(s, d, k, amp, bg, __max_exp_coef);
        return *this;
    }

    ///Read a source with an unspecified PSF.
    SourceInput &SourceInput::operator>>(Core::SourceLocation &source)
    {
        double s, d, k, amp, bg;
        read_source(source, s, d, k, amp, bg);
        return *this;
    }

    ///Read a source with an Elliptical Gaussian PSF that can be 
    ///enabled/disabled according to a column in the input file
    SourceInput &SourceInput::operator>>(IO::OutputSDKSource &source)
    {
        double s, d, k, amp = 1.0, bg = 0.0;
        if(read_source(source, s, d, k, amp, bg)) source.enable();
        else source.disable();
        source.set_psf(s, d, k, amp, bg, __max_exp_coef);
        return *this;
    }

} //End SubPixPhot namespace.
