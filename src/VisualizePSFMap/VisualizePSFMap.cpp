#include "VisualizePSFMap.h"

///Describes the available command line options.
void VisualizePSFMapConfig::describe_options()
{
	_hidden.add_options()
		(
		 "io.psfmap-file",
		 opt::value<std::string>(),
		 "HDF5 file containing the PSF map to visualize."
		);
	_positional.add("io.psfmap-file", 1);

	opt::options_description io_options("Input/Output options");
	io_options.add_options()
		(
		 "io.psfmap-image,o",
		 opt::value<std::string>(),
		 "The name of the output FITS file to create the visualization in."
		)
		(
		 "io.pixels-ds9,p",
		 opt::value<std::string>()->default_value(""),
		 "If passed, the named file is created and filled with ds9 regions "
		 "marking the pixel boundaries in the displayed PSFs."
		)
		(
		 "io.grid-ds9,g",
		 opt::value<std::string>()->default_value(""),
		 "If passed, the named file is created and filled with ds9 regions "
		 "marking the grid boundaries in the displayed PSFs."
		);

	opt::options_description visualization_options(
			"Options defining how the visualization should be performed."
	);
	visualization_options.add_options()
		(
		 "visualize.eval-x",
         opt::value<RealList>(),
		 "The x locations where to evaluate the PSF map (min, max, step)."
		)
		(
		 "visualize.eval-y",
         opt::value<RealList>(),
		 "The y locations where to evaluate the PSF map (min, max, step)."
		)
		(
		 "visualize.psf-x-resolution",
		 opt::value<unsigned>(),
		 "The resolution in the x direction of an individual PSF displayed."
		)
		(
		 "visualize.psf-y-resolution",
		 opt::value<unsigned>(),
		 "The resolution in the y direction to use when displaying the PSF."
		)
		(
		 "visualize.flux",
		 opt::value<double>(),
		 "The flux to give the displayed PSFs."
		)
		(
		 "visualize.bg",
		 opt::value<double>(),
		 "The background to add under the displayed PSFs."
		)
		(
		 "visualize.sdk.min-x",
		 opt::value<double>(),
		 "The minimum offset in x from the center of the PSF to display."
		)
		(
		 "visualize.sdk.max-x",
		 opt::value<double>(),
		 "The maximum offset in x from the center of the PSF to display."
		)
		(
		 "visualize.sdk.min-y",
		 opt::value<double>(),
		 "The minimum offset in y from the center of the PSF to display."
		)
		(
		 "visualize.sdk.max-y",
		 opt::value<double>(),
		 "The maximum offset in y from the center of the PSF to display."
		)
        (
         "visualize.psfmap-variables",
         opt::value<StringList>(),
         "A comma separated list of <variable>=<value> entries defining all "
         "terms in the PSF map which are not x and y."
        );
	_cmdline_config.add(io_options)
				   .add(visualization_options);
}

///Add a single high resolution PSF to the visualization matrix.
template<class EIGEN_MATRIX>
void add_psf(
		///The PSF to visualize
		const PSF &psf,

		///The visualization matrix to add the PSF to.
		EIGEN_MATRIX &visualization,

		///The left boundary of the visualized PSF in the visualization
		///matrix.
		unsigned xmin,

		///The right boundary of the visualized PSF in the visualization
		///matrix.
		unsigned xmax,

		///The bottom boundary of the visualized PSF in the visualizaton
		///matrix.
		unsigned ymin,

		///The top boundary of the visualized PSF in the visualization
		///matrix.
		unsigned ymax,

		///The factor by which to zoom in the PSF in the x direction.
		double xzoom,

		///The factor by which to zoom in the PSF in the x direction.
		double yzoom,

		///The x coordinate of the point relative to which the PSF is
		///defined in the visualization matrix. If NaN, the middle of the x
		///range is used.
		double x0=NaN,

		///The y coordinate of the point relative to which the PSF is
		///defined in the visualization matrix. If NaN, the middle of the y
		///range is used.
		double y0=NaN)
{
	if(std::isnan(x0)) x0 = (xmax + xmin) / 2.0;
	if(std::isnan(y0)) y0 = (ymax + ymin) / 2.0;
	double dxpsf = 1.0 / xzoom, dypsf = 1.0 / yzoom;
	for(unsigned yvis = ymin; yvis < ymax; ++yvis) {
		double ypsf_center = (yvis - y0 + 0.5) / yzoom;
		for(unsigned xvis = xmin; xvis < xmax; ++xvis) {
			double xpsf_center = (xvis - x0 + 0.5) / xzoom;
			visualization(xvis, yvis) = psf.integrate(xpsf_center,
                                                      ypsf_center,
                                                      dxpsf,
                                                      dypsf);
		}
	}
}

///Generates a FITS image of zoomed-in PSFs showing the spatial dependence.
void visualize_map(
		///The map giving how the PSF varies over the image.
		const PSFMap &psf_map,

        ///The arrays of PSF term values to use for evaluating the PSF.
        const std::vector<TermValarray> &expansion_term_arrays,

        ///How many columns to split the displayed PSFs in.
        unsigned num_display_cols,

		///The resolution of each individual diplayed PSF along x.
		unsigned x_resolution,

		///The resolution of each individual diplayed PSF along y.
		unsigned y_resolution,

		///The flux to give each source image.
		double flux,

		///The background to place under the sources.
		double background,

		///The minimum offset in x from the center of the PSF to display.
		double psf_xmin,

		///The maximum offset in x from the center of the PSF to display.
		double psf_xmax,

		///The minimum offset in y from the center of the PSF to display.
		double psf_ymin,

		///The maximum offset in y from the center of the PSF to display.
		double psf_ymax,

		///The name of the file to generate (overwritten if already exists).
		const std::string &outfname)
{
    unsigned num_display_psfs = expansion_term_arrays[0].size(),
             num_display_rows =  num_display_psfs / num_display_cols;
    assert(num_display_psfs % num_display_cols == 0);
    DoubleImageMatrix visualization(num_display_cols * x_resolution,
                                    num_display_rows * y_resolution);

	double xzoom = static_cast<double>(x_resolution) / (psf_xmax - psf_xmin),
		   yzoom = static_cast<double>(y_resolution) / (psf_ymax - psf_ymin),
		   x0_offset = -psf_xmin * xzoom,
           y0_offset = -psf_ymin * yzoom;

    Eigen::VectorXd expansion_terms(expansion_term_arrays.size());
    for(unsigned psf_index = 0; psf_index < num_display_psfs; ++psf_index) {

        for(unsigned term_i = 0; term_i < expansion_terms.size(); ++term_i)
            expansion_terms[term_i] =
                expansion_term_arrays[term_i][psf_index];
        PSF *psf=psf_map(expansion_terms);

        unsigned display_row = psf_index / num_display_cols,
                 display_col = psf_index % num_display_cols;
        double x_min = display_col * x_resolution,
               y_min = display_row * y_resolution;
        add_psf(*psf,
                visualization,
                x_min, x_min + x_resolution,
                y_min, y_min + y_resolution,
                xzoom, yzoom,
                x0_offset + x_min, y0_offset + y_min);
        delete psf;
    }
	visualization *= flux / xzoom / yzoom;
	visualization.array() += background / xzoom / yzoom;
	save_matrix_as_fits(outfname, visualization);
}

///Saves a given list of boundaries to a file.
void save_boundaries(
		///The coordinates of the vertical boundaries to display in PSF
		///coordinates.
		const std::vector<double> x_boundaries,

		///The coordinates of the horizontal boundaries to display in PSF
		///coordinates.
		const std::vector<double> y_boundaries,

		///The number of different x positions to sample the PSF at.
		unsigned x_splits,

		///The number of different x positions to sample the PSF at.
		unsigned y_splits,

		///The resolution of each individual diplayed PSF along x.
		unsigned x_resolution,

		///The resolution of each individual diplayed PSF along y.
		unsigned y_resolution,

		///The minimum offset in x from the center of the PSF displayed.
		double psf_xmin,

		///The maximum offset in x from the center of the PSF displayed.
		double psf_xmax,

		///The minimum offset in y from the center of the PSF displayed.
		double psf_ymin,

		///The maximum offset in y from the center of the PSF displayed.
		double psf_ymax,

		///The color to use for the regions.
		const std::string &color,

		///The name of the file to create/overwrite.
		const std::string &filename)
{
	std::ofstream outf(filename.c_str());
	outf << "global color="+color+" dashlist=3 8 width=1 font=\"helvetica 10"
		" normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1"
		" include=1 source=1" << std::endl;
	int xmax=x_splits*x_resolution, ymax=y_splits*y_resolution;
	double x_scale=x_resolution/(psf_xmax-psf_xmin),
		   y_scale=y_resolution/(psf_ymax-psf_ymin);
	for(unsigned x_ind=0; x_ind<x_splits; ++x_ind)
		for(std::vector<double>::const_iterator
				xb_i=x_boundaries.begin();
				xb_i!=x_boundaries.end(); ++xb_i) {
			double image_x=x_ind*x_resolution
						   +
						   x_scale*(*xb_i-psf_xmin)
;
			outf << "physical; line(" << image_x << ",1," << image_x << ","
				<< ymax << ")# line=0 0" << std::endl;
		}
	for(unsigned y_ind=0; y_ind<y_splits; ++y_ind)
		for(std::vector<double>::const_iterator
				yb_i=y_boundaries.begin();
				yb_i!=y_boundaries.end(); ++yb_i) {
			double image_y=y_ind*y_resolution
						   +
						   y_scale*(*yb_i-psf_ymin);
			outf << "physical; line(1," << image_y << "," << xmax << ","
				<< image_y << ")# line=0 0" << std::endl;
		}
	outf.close();
}

void create_pixel_regions(
		///The number of different x positions to sample the PSF at.
		unsigned x_splits,

		///The number of different x positions to sample the PSF at.
		unsigned y_splits,

		///The resolution of each individual diplayed PSF along x.
		unsigned x_resolution,

		///The resolution of each individual diplayed PSF along y.
		unsigned y_resolution,

		///The minimum offset in x from the center of the PSF displayed.
		double psf_xmin,

		///The maximum offset in x from the center of the PSF displayed.
		double psf_xmax,

		///The minimum offset in y from the center of the PSF displayed.
		double psf_ymin,

		///The maximum offset in y from the center of the PSF displayed.
		double psf_ymax,

		///The name of the file to create/overwrite.
		const std::string &filename)
{
	int xmin=static_cast<int>(std::ceil(psf_xmin)),
		xmax=static_cast<int>(std::floor(psf_xmax)),
		ymin=static_cast<int>(std::ceil(psf_ymin)),
		ymax=static_cast<int>(std::floor(psf_ymax));
	std::vector<double> x_boundaries(xmax-xmin+1), y_boundaries(ymax-ymin+1);
	for(int i=0; i<static_cast<int>(x_boundaries.size()); ++i)
		x_boundaries[i]=xmin+i;
	for(int i=0; i<static_cast<int>(y_boundaries.size()); ++i)
		y_boundaries[i]=ymin+i;
	save_boundaries(x_boundaries,
					y_boundaries,
					x_splits,
					y_splits,
					x_resolution,
					y_resolution,
					psf_xmin,
					psf_xmax,
					psf_ymin,
					psf_ymax,
					"green",
					filename);
}

void create_grid_regions(
		///The PSF map being visualized.
		const PiecewiseBicubicPSFMap &psf_map,

		///The number of different x positions to sample the PSF at.
		unsigned x_splits,

		///The number of different x positions to sample the PSF at.
		unsigned y_splits,

		///The resolution of each individual diplayed PSF along x.
		unsigned x_resolution,

		///The resolution of each individual diplayed PSF along y.
		unsigned y_resolution,

		///The name of the file to create/overwrite.
		const std::string &filename)
{
	const std::vector<double> &x_grid=psf_map.x_grid(),
							  &y_grid=psf_map.y_grid();
	save_boundaries(x_grid,
					y_grid,
					x_splits,
					y_splits,
					x_resolution,
					y_resolution,
					x_grid.front(),
					x_grid.back(),
					y_grid.front(),
					y_grid.back(),
					"red",
					filename);
}

///Return the name and value in an expression like \<name\>=\<value\>
PSFMapVariableType parse_variable_expression(
        const std::string &expression,
        unsigned array_size
)
{
    size_t equal_pos = expression.find_first_of('=');
    double value;
    std::istringstream(expression.substr(equal_pos + 1)) >> value;
    return PSFMapVariableType(expression.substr(0, equal_pos),
                              std::valarray<double>(value, array_size));
}

///Fill the variables of the PSF map and return how many columns to display.
///
///The values are ordered as they should be displayed in the map: with the
///lower left cornern first with each row's entries sequential (row-major
///format).
unsigned fill_psfmap_variables(
        ///The command line options.
        const VisualizePSFMapConfig &options,

        ///The object to fill with the variables.
        PSFMapVarListType &variables
)
{
    RealList::const_iterator
        eval_iter = options["visualize.eval-x"].as<RealList>().begin();
    double x_min = *eval_iter++,
           x_max = *eval_iter++,
           x_step = *eval_iter++;
    assert(eval_iter == options["visualize.eval-x"].as<RealList>().end());
    eval_iter = options["visualize.eval-y"].as<RealList>().begin();
    double y_min = *eval_iter++,
           y_max = *eval_iter++,
           y_step = *eval_iter++;
    assert(eval_iter == options["visualize.eval-y"].as<RealList>().end());
    assert(x_max > x_min);
    assert(y_max > y_min);

    unsigned num_display_cols = unsigned((x_max - x_min) / x_step) + 1,
             num_display_psfs = num_display_cols
                                *
                                (unsigned((y_max - y_min) / y_step) + 1);
    variables.push_back(
            PSFMapVariableType("x", std::valarray<double>(num_display_psfs))
    );
    std::valarray<double>& eval_x = variables.back().second;
    variables.push_back(
            PSFMapVariableType("y", std::valarray<double>(num_display_psfs))
    );
    std::valarray<double>& eval_y = variables.back().second;
    unsigned display_index = 0;
    for(double y = y_min; y < y_max; y += y_step)
        for(double x = x_min; x < x_max; x += x_step) {
            assert(display_index < num_display_psfs);
            eval_x[display_index] = x;
            eval_y[display_index] = y;
            ++display_index;
        }

    const StringList &psf_var_strings =
        options["visualize.psfmap-variables"].as<StringList>();
    for(
            std::list<std::string>::const_iterator
                var_str_i = psf_var_strings.begin();
            var_str_i != psf_var_strings.end();
            ++var_str_i
    )
        variables.push_back(parse_variable_expression(*var_str_i,
                                                      num_display_psfs));
    return num_display_cols;
}

///\brief Fill the terms on which the PSF map depends and return how many
///columns to display.
///
///The values are ordered as they should be displayed in the map: with the
///lower left cornern first with each row's entries sequential (row-major
///format).
unsigned fill_psfmap_terms(
        ///The command line options.
        const VisualizePSFMapConfig &options,

        ///The expression defining the terms which participate in the PSF
        ///map.
        const std::string &expansion_term_expression,

        ///The object to fill with PSF terms.
        std::vector<TermValarray> &expansion_term_arrays
)
{
    PSFMapVarListType psfmap_variables;
    unsigned num_display_cols = fill_psfmap_variables(options,
                                                      psfmap_variables);
    evaluate_term_expression(expansion_term_expression,
                             psfmap_variables.begin(),
                             psfmap_variables.end(),
                             expansion_term_arrays);
    return num_display_cols;
}

int main(int argc, char *argv[])
{
#ifndef DEBUG
	try {
#endif
		VisualizePSFMapConfig options(argc, argv);
		if(!options.proceed()) return 1;
		H5IODataTree psfmap_data;
		std::set<std::string> required_quantities(
				PiecewiseBicubicPSFMap::required_data_tree_quantities()
		);
		required_quantities.insert(
				EllipticalGaussianPSFMap::required_data_tree_quantities()
																.begin(),
				EllipticalGaussianPSFMap::required_data_tree_quantities()
																.end()
		);
		psfmap_file.read(required_quantities.begin(),
						 required_quantities.end(),
						 psfmap_data,
						 false);
		std::string psf_model=psfmap_data.get<std::string>("psffit.model",
														   "",
														   translate_string);
        std::vector<TermValarray> expansion_term_arrays;
        unsigned
            num_display_cols = fill_psfmap_terms(
                    options,
                    psfmap_data.get<std::string>("psffit.terms",
                                                 "",
                                                 translate_string),
                    expansion_term_arrays
            ),
            num_display_rows = (expansion_term_arrays[0].size()
                                /
                                num_display_cols);
        unsigned psf_x_resolution = options[
                     "visualize.psf-x-resolution"
                 ].as<unsigned>(),
                 psf_y_resolution = options[
                     "visualize.psf-y-resolution"
                 ].as<unsigned>();

		PSFMap *psf_map;
		double psf_xmin, psf_xmax, psf_ymin, psf_ymax;
		if(psf_model=="bicubic") {
			PiecewiseBicubicPSFMap
				*bicubic_map=new PiecewiseBicubicPSFMap(psfmap_data);
			psf_map=bicubic_map;
			psf_xmin=bicubic_map->x_grid().front();
			psf_xmax=bicubic_map->x_grid().back();
			psf_ymin=bicubic_map->y_grid().front();
			psf_ymax=bicubic_map->y_grid().back();
			std::string
				grid_regions_fname=options["io.grid-ds9"].as<std::string>();
			if(grid_regions_fname != "")
                create_grid_regions(*bicubic_map,
                                    num_display_cols, num_display_rows,
                                    psf_x_resolution, psf_y_resolution,
                                    grid_regions_fname);
		} else  {
			assert(psf_model=="sdk");
			psf_map=new EllipticalGaussianPSFMap(psfmap_data);
			psf_xmin=options["visualize.sdk.min-x"].as<double>();
			psf_xmax=options["visualize.sdk.max-x"].as<double>();
			psf_ymin=options["visualize.sdk.min-y"].as<double>();
			psf_ymax=options["visualize.sdk.max-y"].as<double>();
		}
		visualize_map(*psf_map,
                      expansion_term_arrays,
                      num_display_cols,
                      psf_x_resolution, psf_y_resolution,
                      options["visualize.flux"].as<double>(),
                      options["visualize.bg"].as<double>(),
                      psf_xmin, psf_xmax, psf_ymin, psf_ymax,
                      options["io.psfmap-image"].as<std::string>());
		std::string
			pixel_regions_fname=options["io.pixels-ds9"].as<std::string>();
		if(pixel_regions_fname != "")
            create_pixel_regions(num_display_cols, num_display_rows,
                                 psf_x_resolution, psf_y_resolution,
                                 psf_xmin, psf_xmax, psf_ymin, psf_ymax,
                                 pixel_regions_fname);
		delete psf_map;
		return 0;
#ifndef DEBUG
	} catch(Error::General &ex) {
		std::cerr << ex.what() << ":" << ex.get_message() << std::endl;
		return 1;
	}
#endif
}
