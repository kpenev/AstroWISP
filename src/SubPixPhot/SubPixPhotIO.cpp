#include "SubPixPhotIO.h"

namespace SubPixPhot {

    static void output_double(double          value,
                              std::ostream&   os,
                              int             precision   = 3,
                              int             width       = 10
    )
    {
        os << std::setfill(' ')
           << std::setprecision(precision)
           << std::fixed
           << std::setw(width)
           << value;
    }

    static void output_int(int value, std::ostream& os, int width = 5)
    {
        os << std::setfill(' ') << std::setw(width) << value;
    }

    void output_to_stdout(
        const std::list<IO::OutputSDKSource>&   sources,
        const std::list<Phot::Columns>&         columns,
        double                                  mag_1ADU,
        double                                  gain,
        std::ostream&                           os /*= std::cout*/
    )
    {
        unsigned flux_ind = 0,
                 flux_err_ind = 0,
                 mag_ind = 0,
                 mag_err_ind = 0,
                 flag_ind = 0;
        os << '#';
        for(
            std::list<Phot::Columns>::const_iterator col_i = columns.begin();
            col_i != columns.end();
            ++col_i
        ) {
            int corr = (col_i == columns.begin() ? 1 : 0);
            switch(*col_i) {
                case Phot::id : os << std::setw(15-corr) << "ID"; break;
                case Phot::x : os << std::setw(10-corr) << "x"; break;
                case Phot::y : os << std::setw(10-corr) << "y"; break;
                case Phot::S : os << std::setw(10-corr) << "S"; break;
                case Phot::D : os << std::setw(10-corr) << "D"; break;
                case Phot::K : os << std::setw(10-corr) << "K"; break;
                case Phot::A : os << std::setw(18-corr) << "Amp"; break;
                case Phot::bg : os << std::setw(15-corr) << "Bg"; break;
                case Phot::bg_err : os << std::setw(15-corr)
                                    << "BgErr"; break;
                case Phot::chi2 : os << std::setw(15-corr) << "Chi^2"; break;
                case Phot::sn : os << std::setw(15-corr) << "S/N"; break;
                case Phot::npix : os << std::setw(5-corr) << "Npix"; break;
                case Phot::nbgpix: os << std::setw( 5 - corr ) << "BgPix"; 
                                   break;
#ifdef DEBUG
                case Phot::time: os << std::setw(15-corr) << "time"; break;
#endif
                case Phot::flux : os << std::setw(16-corr) << "Flux["
                                     << std::setfill('0')
                                     << std::setw(3) << flux_ind++ << "]"
                                     << std::setfill(' ');
                                  break;
                case Phot::flux_err : os << std::setw(16) << "Flux error["
                                         << std::setfill('0')
                                         << std::setw(3) << flux_err_ind++
                                         << ']'
                                         << std::setfill(' ');
                                      break;
                case Phot::mag : os << std::setw(6-corr) << "Mag["
                                    << std::setfill('0')
                                    << std::setw(3) << mag_ind++
                                    << ']'
                                    << std::setfill(' ');
                                 break;
                case Phot::mag_err : os << std::setw(8-corr) << "MagErr["
                                        << std::setfill('0')
                                        << std::setw(3) << mag_err_ind++
                                        << ']'
                                        << std::setfill(' ');
                                     break;
                case Phot::flag : os << std::setw(5-corr) << "Flag["
                                     << std::setfill('0')
                                     << std::setw(3) << flag_ind++
                                     << ']'
                                     << std::setfill(' ');
                                  break;
                default : os << "?????";
            }
            os << ' ';
        }
        os << std::endl;
        for(
            std::list<IO::OutputSDKSource>::const_iterator 
                src_i = sources.begin();
            src_i != sources.end();
            ++src_i
        ) {
            flux_ind = 0;
            flux_err_ind = 0;
            mag_ind = 0;
            mag_err_ind = 0;
            flag_ind = 0;
            for(
                std::list<Phot::Columns>::const_iterator
                    col_i = columns.begin();
                col_i != columns.end();
                ++col_i
            ) {
                switch(*col_i) {
                    case Phot::id: os << src_i->id(); break;
                    case Phot::x: output_double(src_i->x(), os); break;
                    case Phot::y: output_double(src_i->y(), os); break;
                    case Phot::S: output_double(src_i->psf_s(), os); break;
                    case Phot::D: output_double(src_i->psf_d(), os); break;
                    case Phot::K: output_double(src_i->psf_k(), os); break;
                    case Phot::A: output_double(
                                       src_i->psf_amplitude() / gain, 
                                       os,
                                       3,
                                       18
                                   );
                                   break;
                    case Phot::bg: output_double(
                                       src_i->background().value()/gain,
                                       os,
                                       3,
                                       15
                                   );
                                   break;
                    case Phot::bg_err: output_double(
                                           src_i->background().error() / gain,
                                           os,
                                           3,
                                           15
                                       );
                                       break;
                    case Phot::chi2: output_double(src_i->reduced_chi2(),
                                                   os,
                                                   3,
                                                   15);
                                     break;
                    case Phot::sn: output_double(src_i->signal_to_noise(),
                                                 os,
                                                 2,
                                                 15);
                                   break;
                    case Phot::npix: output_int(src_i->pixel_count(), os);
                                     break;
                    case Phot::nbgpix: output_int(
                                           src_i->background().pixels(),
                                           os
                                       );
                                       break;
#ifdef DEBUG
                    case Phot::time: output_double(src_i->processing_time(),
                                                   os,
                                                   3,
                                                   15);
                                     break;
#endif
                    case Phot::flux : output_double(
                                          src_i->flux()[flux_ind++].value()
                                          /
                                          gain, 
                                          os,
                                          3,
                                          20
                                      );
                                      break;
                    case Phot::flux_err :
                        output_double(
                            src_i->flux()[flux_err_ind++].error() / gain,
                            os,
                            3,
                            20
                        );
                        break;
                    case Phot::mag:
                    {
                        const double mag = (
                            mag_1ADU
                            -
                            2.5 * log10(
                                src_i->flux()[mag_ind++].value() / gain
                            )
                        );
                        output_double(mag, os);
                    }
                    break;
                    case Phot::mag_err : 
                    {
                        output_double(
                            -2.5 * log10(1.0
                                         -
                                         src_i->flux()[mag_err_ind].error()
                                         /
                                         src_i->flux()[mag_err_ind].value()),
                            os,
                            3, 
                            12
                        );
                        ++mag_err_ind;
                    }
                    break;
                    case Phot::flag: 
                        switch(src_i->flux()[flag_ind++].flag()) {
                            case Core::GOOD : os << "    G    "; break;
                            case Core::SATURATED : os << "    C    "; break;
                            case Core::BAD : os << "    X    "; break;
                            default : os << "    ?    "; break;
                        }
                        break;
                    default : os << "?????";
                }
                os << ' ';
            }
            os << std::endl;
        }
    }

    void write_header(std::ostream      &os,
                      int               argc,
                      char**            argv,
                      const std::string &version_string)
    {
#ifdef HATP_SVN
        os << "# HATpipe SVN revision: " << HATP_SVN << std::endl;
#endif
        if(version_string.size()>0)
            os << "# Created by " << version_string << std::endl;
        if(argc>0) {
            os << "# Command line:";
            for(int i=0; i<argc; i++) os << " '" << argv[i] << "'";
            os << std::endl;
        }
    }

    std::list<double> read_sdk_coef(const std::string &fname)
    {
        std::list<double> result;
        if(fname == "") return result;
        std::ifstream infile(fname.c_str());
        char c;
        std::string str;
        infile >> c;
        assert(c == '#');
        infile >> str;
        assert(std::string(str) == "PSF");
        infile >> str;
        assert(std::string(str) == "coefficients:");
        std::getline(infile, str);
        for(std::istringstream line(str); !line.eof();) {
            double coef;
            line >> coef;
            if(!line) throw Error::IO(
                "Failed to parse initial guess coefficients from '"
                +
                str
                +
                "'."
            );
            result.push_back(coef);
        }
        return result;
    }

} //End SubPixPhot namespace.
