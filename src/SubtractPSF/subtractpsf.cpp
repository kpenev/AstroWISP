
#include <iostream>

#include "subtractpsf.h"


void SubtractPSFCfg::describe_options()
{
	_hidden.add_options()
		(
		 "io.image", 
		 opt::value<std::string>(),
		 "Reduced file to fit the PSF of."
		)
		(
		 "io.psffit",
		 opt::value<std::string>(),
		 "The file produced by FitPSF containing the PSF map and fluxes to "
		 "subtract."
		);
	_positional.add("io.image", 1)
			   .add("io.psffit", 1);
	opt::options_description io_options("Input/Output options");
	io_options.add_options()
		(
		 "io.psfsubtracted-image,o",
		 opt::value<std::string>(),
		 "Filename for the PSF subtracted image."
		)
		(
		 "io.subpix,s",
		 opt::value<std::string>()->default_value(""),
		 "FITS file representing the sensitivity structure of a pixel. "
		 "Uniform sensitivity is assumed if empty string or omitted."
		);
	opt::options_description background_options(
			"Options defining how to measure the background behing the "
			"sources.");
	background_options.add_options()
		(
		 "bg.annulus,b",
		 opt::value<BackgroundAnnulus>(),
		 "Specifies that an annulus with the given inner radius centered "
		 "around the source should be used to estimate the background and "
		 "its error."
		);
	_cmdline_config.add(io_options)
				   .add(background_options);
}

std::string SubtractPSFCfg::usage_help(const std::string &prog_name) const
{
	return "Usage: " + prog_name + " [OPTION ...] <FITS> <HDF5>";
}

// ----------------------------------------------------------------------------
// local definitions

typedef FitsImage<double>   FITSImage;
typedef std::vector<double> dvector;

static void subtract(
	const H5IODataTree&				psf_data,
    FITSImage&                      img,
    const FITSImage*                subpixmapptr,
    const std::string&              output_fname,
	double							max_psf_extent
);

static void subtract_src(
    FITSImage&                      img,
    const PSF*                      psfmap,
    double                          xcord,
    double                          ycord,
    double                          flux,
    double                          xminus,
    double                          xplus,
    double                          yminus,
    double                          yplus,
    const dvector&                  subpixvals,
    int                             x_subgridsize,
    int                             y_subgridsize
);


// ----------------------------------------------------------------------------

int main( int argc, char** argv )
{
    // parse command-line options
    SubtractPSFCfg options( argc, argv );
    if (!options.proceed()) return 1;

    H5::Exception::dontPrint();

    // open the original FITS image
#ifndef DEBUG	
    try {
#endif
#ifdef TRACK_PROGRESS
		std::cerr << "Reading sources to subtract and PSF map." << std::endl;
#endif
        // open HDF5 containing the PSF map
        SubPixHDF5File hdf5file(
				options["io.psffit"].as<std::string>().c_str(),
				H5F_ACC_RDONLY
		);

		H5::Group fitpsf=hdf5file.openGroup("PSFFit");

        // read PSF map data and quantities required to get the amplitudes.
		std::set<std::string> quantities_to_read(
				PiecewiseBicubicPSFMap::required_data_tree_quantities()
		);
		quantities_to_read.insert("psffit.flux");
        quantities_to_read.insert("psffit.tool");
		quantities_to_read.insert("psffit.magnitude_1adu");
		quantities_to_read.insert("psffit.mag");
		quantities_to_read.insert("projsrc.x");
		quantities_to_read.insert("projsrc.y");
		H5IODataTree psf_data;
		hdf5file.read(quantities_to_read.begin(),
					  quantities_to_read.end(),
					  psf_data,
					  false);
        PSF::fill_psf_fluxes(psf_data);

        // input image
        FITSImage img( options["io.image"].as<std::string>() );

        // output image
        FITSImage outimg = img;

        // subpixel map if any
        FITSImage subpixmap;
        FITSImage* subpixmapptr = 0;
        const std::string& 
			subpixfname = options["io.subpix"].as<std::string>();
        if ( !subpixfname.empty() ) {
            subpixmap.open( subpixfname );
            subpixmapptr = &subpixmap;
        }

#ifdef TRACK_PROGRESS
		std::cerr << "Starting subtraction." << std::endl;
#endif
        // do the work (outimg: input/output)
        subtract(
				psf_data,
				outimg,
				subpixmapptr,
				options["io.psfsubtracted-image"].as<std::string>(),
				options["bg.annulus"].as<BackgroundAnnulus>().inner_radius()
		);
#ifdef TRACK_PROGRESS
		std::cerr << "Done." << std::endl;
#endif
#ifndef DEBUG
    }
    catch ( const char* msg ) {
        std::cerr
            << argv[0]
            << ": "
            << msg
            << std::endl;

        return 2;
    }
    catch ( const H5::Exception& ex ) {
        std::cerr
            << argv[0]
            << ": cannot open file "
            << options["io.psffit"].as<std::string>()
            << " ("
            << ex.getCDetailMsg()
            << ")"
            << std::endl;
        ;
    }
    catch ( Error::General& ex ) {
        std::cerr
            << argv[0]
            << ": "
            << ex.what()
            << ":"
            << ex.get_message()
            << std::endl;
    }
    catch ( ... ) {
        std::cerr
            << argv[0]
            << ": could not open FITS file "
            << options["io.image"].as<std::string>()
            << std::endl;
        ;

        return 1;
    }
#endif

    return 0;
}

// ----------------------------------------------------------------------------

void subtract(
	const H5IODataTree&				psf_data,
    FITSImage&                      img,
    const FITSImage*                subpixmapptr,
    const std::string&              output_fname,
	double 							max_psf_extent
)
{
	PSFMap *psfmap;
	std::string psf_model=psf_data.get<std::string>("psffit.model",
													"",
													translate_string);
    std::string fit_tool=psf_data.get<std::string>("psffit.tool",
                                                   "",
                                                   translate_string);
	double xminus, xplus, yminus, yplus;
	if(psf_model=="bicubic") {
		PiecewiseBicubicPSFMap 
			*bicubic_map=new PiecewiseBicubicPSFMap(psf_data);
		const dvector& xgrid = bicubic_map->x_grid();
		const dvector& ygrid = bicubic_map->y_grid();
		xminus = xgrid[0];
		xplus = xgrid.back();
		yminus = ygrid[0];
		yplus = ygrid.back();
		psfmap=bicubic_map;
	} else {
		assert(psf_model=="sdk");
		psfmap=new EllipticalGaussianPSFMap(psf_data);
		xminus = -max_psf_extent;
		xplus = max_psf_extent;
		yminus = -max_psf_extent;
		yplus = max_psf_extent;
	}

    int x_subgridsize = 1;
    int y_subgridsize = 1;
    dvector subpixvals;

    if ( fit_tool == "FitPRF" ) x_subgridsize = y_subgridsize = 0;
    else {
        assert( fit_tool == "FitPSF" );
        if ( subpixmapptr ) {
            x_subgridsize = subpixmapptr->data().x_resolution();
            y_subgridsize = subpixmapptr->data().y_resolution();

            subpixvals.resize( x_subgridsize * y_subgridsize );
            dvector::iterator it = subpixvals.begin();
            for ( int yy = 0; yy < y_subgridsize; ++yy ) {
                for ( int xx = 0; xx < x_subgridsize; ++xx, ++it ) {
                    *it = subpixmapptr->data()( xx, yy );
                }
            }
        } else {
            // no subpixel map
            subpixvals.push_back( 1 );
        }
    }

	OutputArray<double> x(psf_data.get<boost::any>("projsrc.x")),
						y(psf_data.get<boost::any>("projsrc.y")),
						flux_values(psf_data.get<boost::any>("psffit.flux"));
	size_t num_sources=x.size();

	assert(y.size()==num_sources);
	assert(flux_values.size()==num_sources);
    for (size_t source_index=0; source_index<num_sources; ++source_index) {
        double xcord = x[source_index];
        double ycord = y[source_index];
        double flux = flux_values[source_index];

        PSF* psf = (*psfmap)( xcord, ycord );

#ifdef TRACK_PROGRESS
		std::cerr << "Subtracting source " << source_index << "/"
				  << num_sources << ": x=" << xcord << ", y=" << ycord
				  << ", flux=" << flux << std::endl;
#endif

        subtract_src(
            img,
            psf,
            xcord,
            ycord,
            flux,
            xminus,
            xplus,
            yminus,
            yplus,
            subpixvals,
            x_subgridsize,
            y_subgridsize
        );

        delete psf;
    }
#ifdef TRACK_PROGRESS
		std::cerr << "Saving subtracted image." << std::endl;
#endif
    img.save( output_fname );
}

// ----------------------------------------------------------------------------

void subtract_src(
    FITSImage&      img,
    const PSF*      psf,
    double          xcord,
    double          ycord,
    double          flux,
    double          xminus,
    double          xplus,
    double          yminus,
    double          yplus,
    const dvector&  subpixvals,
    int             x_subgridsize,
    int             y_subgridsize
)
{
    if ( x_subgridsize == 0 ) {
        assert( y_subgridsize == 0 );
        x_subgridsize = 1;
        y_subgridsize = 1;
        xminus += 0.5;
        yminus += 0.5;
        xplus -= 0.5;
        yplus -= 0.5;
    }
    int naxis1 = img.data().x_resolution() - 1; // last valid index
    int naxis2 = img.data().y_resolution() - 1; // last valid index

    int xmin = std::min( naxis1, std::max( 0, int( xcord + xminus ) ) );
    int xmax = std::max( 0, std::min( naxis1, int( xcord + xplus ) ) );
    int ymin = std::min( naxis2, std::max( 0, int( ycord + yminus ) ) );
    int ymax = std::max( 0, std::min( naxis2, int( ycord + yplus ) ) );

    int xsize = xmax - xmin + 1;
    int ysize = ymax - ymin + 1;

    double xgridsize = 1.0 / x_subgridsize;
    double ygridsize = 1.0 / y_subgridsize;

    dvector psfvals( xsize * ysize, 0 );
    double psfsum = 0;

    for ( int yy = ymin; yy <= ymax; ++yy ) {
        int valoffset = ( yy - ymin ) * xsize;

        for ( int ystep = 0; ystep < y_subgridsize; ++ystep ) {
            double yoff = yy + ( ystep + 0.5 ) * ygridsize;

            int subpixoffset = ystep * x_subgridsize;

            for ( int xx = xmin; xx <= xmax; ++xx ) {
                for ( int xstep = 0; xstep < x_subgridsize; ++xstep ) {
                    double xoff = xx + ( xstep + 0.5 ) * xgridsize;
                    double psfval, subpixval;
                    if ( subpixvals.size() == 0 ) {
                        psfval = (*psf)( xoff - xcord, yoff - ycord );
                        subpixval = 1;
                    } else {
                        psfval = psf->integrate(
                            xoff - xcord,
                            yoff - ycord,
                            xgridsize,
                            ygridsize
                        );

                        subpixval = subpixvals[subpixoffset + xstep];
                    }

                    psfvals[valoffset + xx - xmin] += psfval * subpixval;
                    psfsum += psfval;
                }
            }
        }
    }

    dvector::iterator it = psfvals.begin();
    for ( int yy = ymin; yy <= ymax; ++yy ) {
        for ( int xx = xmin; xx <= xmax; ++xx, ++it ) {
            *it *= flux / psfsum;

            img.data()( xx, yy ) -= *it;
        }
    }
}

// ----------------------------------------------------------------------------
