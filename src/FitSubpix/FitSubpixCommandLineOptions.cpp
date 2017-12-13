#include "FitSubpixCommandLineOptions.h"
#include <argtable2.h>
#include <sstream>
#include <gsl/gsl_multimin.h>
#include <iomanip>

FitSubpixCommandLineOptions::FitSubpixCommandLineOptions(int argc,
		char **argv)
{
	__parsed_ok=false;
	arg_file *input_files_option = arg_filen(NULL, NULL, 
			"<filename>", 2, 10000, 
			"A list of fits filenames followed by the exact same number of "
			"source filenames to be used for photometry (the variance of "
			"which will be minimized). All source list files should contain "
			"exactly the same source IDs in the same order. No more than "
			"5000 fits files can be used for fitting. WARNING: All input "
			"frames and sources are held in memory.");
	arg_int *x_split_option = arg_int0("x", "x-split", "<int>",
			"The number of pieces in which to split pixels in the x "
			"direction. Default: 4.");
	arg_int *y_split_option = arg_int0("y", "y-split", "<int>",
			"The number of pieces in which to split pixels in the y "
			"direction. Default: 4.");
	arg_dbl *aperture_option = arg_dbl0("a", "aperture", "<real value>",
			"The aperture to used to do photometry. Default: 3.0");
	arg_str *bgannulus_option	= arg_str0("b", "bg-annulus", 
			"<inner radius>,<width>", "Specifies that an annulus with the "
			"given inner radius centered around the source should be used to"
			" estimate the background and its error. Default: 6,7");
    arg_str* incolumns_option = arg_str0(
        NULL,
        "input-columns", "<column>,[<column>,[...]]",
        "A comma separated list of input column names. Recognized values: "
        "id, x, y, S, D, K, A|amp, flux, flux_err, mag, mag_err, flag, bg, "
        "bg_err, sn, npix, nbgpix. Where a group of consecutive columns that "
        "depend on aperture are repeated automatically. Default: "
        "'id, x, y, S, D, K, A|amp."
    );
	arg_str *fit_method_option = arg_str0("m", "fit-method", 
			"<simplex|multinest|nr|siman>",
			"The method to use to fit for the subpixel structure. By default"
			" uses Newton-Raphson (nr).");
	arg_dbl *siman_boltzman_k = arg_dbl0(NULL, "siman-boltzman",
			"<dobule>", 
			"The boltzman constant (default 1e7) to use for the simulated "
			"annealing fitting method.");
	arg_dbl *siman_initial_temperature = arg_dbl0(NULL,
			"siman-start-temperature", "<double>", 
			"The initial temperature (default 1.0) to use for the simulated "
			"annealing fitting method.");
	arg_dbl *siman_cooling_rate = arg_dbl0(NULL, "siman-cooling-rate",
			"<double>",	"The cooling rate (default 1.00001) to use for the "
			"simulated annealing fitting method.");
	arg_dbl *siman_min_temperature = arg_dbl0(NULL, "siman-min-temperature",
			"<double>", "The temperature at which simulated annealing "
			"stops. Default: 1e-3.");
	arg_int *siman_ntries = arg_int0(NULL, "siman-ntries", "<unsigned int>",
			"The number of points to try for each simulated annealing step."
			" Default: 100.");
	arg_dbl *siman_max_step = arg_dbl0(NULL, "siman-max-step", "<double>",
			"The maximum simulated annealing step. Default 1.");
	arg_dbl *max_exp_coef = arg_dbl0(NULL, "max-exp-coef", "<double>",
			"The maximum allowed value for any term appearing in an exponent"
			"while integrating PSFs. Larger values generally result in "
			"faster code, but could lead to extreme numerical round-off "
			"error.");
	arg_lit *help_option = arg_lit0("h", "help", 
					"Print this help and exit.");
	struct arg_end *end = arg_end(100);
	void *argtable[] = {input_files_option, x_split_option, y_split_option, 
						aperture_option, bgannulus_option, incolumns_option,
						fit_method_option, siman_boltzman_k,
						siman_initial_temperature, siman_cooling_rate, 
						siman_min_temperature, siman_ntries, siman_max_step,
						max_exp_coef, help_option, end};
	if(arg_nullcheck(argtable) != 0) 
		throw Error::CommandLine("Failed to allocate argument table.");

	x_split_option->ival[0]			= 4;
	y_split_option->ival[0]			= 4;
	aperture_option->dval[0]		= 3.0;
	bgannulus_option->sval[0]		= "6,7";
	incolumns_option->sval[0]		= "id,x,y,S,D,K,A";
	fit_method_option->sval[0]		= "nr";
	siman_boltzman_k->dval[0]		= 1e7;
	siman_initial_temperature->dval[0] = 1.0;
	siman_cooling_rate->dval[0]		= 1.00001;
	siman_min_temperature->dval[0]	= 1e-3;
	siman_ntries->ival[0]			= 100;
	siman_max_step->dval[0]			= 1.0;
	max_exp_coef->dval[0]			= 1;

	int nerrors=arg_parse(argc, argv, argtable);
	if(input_files_option->count%2) nerrors++;
	if(help_option->count>0 || nerrors>0) {
        printf("Usage: %s", "FitSubpix");
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		if(help_option->count==0) arg_print_errors(stdout, end, "SubPixPhot");
		if(input_files_option->count%2)
			printf("Odd number of positional arguments found!\n");
		return;
	}
	__x_split=x_split_option->ival[0];
	__y_split=y_split_option->ival[0];
	__aperture=aperture_option->dval[0];
	std::list<double> bgan_values=parse_real_list(bgannulus_option->sval[0], 
					"--bg-annulus", 2, 2);
	__bg_annulus.inner_radius()=bgan_values.front();
	__bg_annulus.width()=bgan_values.back();
	int nframes=input_files_option->count/2;
	__frames.resize(nframes);
	__sources.resize(nframes);
	for(int i=0; i<nframes; i++) {
		__frames[i]=input_files_option->filename[i];
		__sources[i]=input_files_option->filename[nframes+i];
	}
	__input_columns=parse_column_list(incolumns_option->sval[0], 
					aperture_option->count, "--input-columns");

	__siman_ntries=siman_ntries->ival[0];

	__siman_max_step=siman_max_step->dval[0];

	__siman_boltzman_k=siman_boltzman_k->dval[0];

	__siman_start_temperature=siman_initial_temperature->dval[0];

	__siman_cooling_rate=siman_cooling_rate->dval[0];

	__siman_min_temperature=siman_min_temperature->dval[0];

	__max_exp_coef=max_exp_coef->dval[0];

	verify_input_columns(__input_columns);
	std::string fit_method(fit_method_option->sval[0]);
	if(fit_method=="simplex") __fit_method=GSL_simplex;
	else if(fit_method=="multinest") __fit_method=MultiNest;
	else if(fit_method=="nr") __fit_method=NewtonRaphson;
	else if(fit_method=="siman") __fit_method=GSL_simulated_annealing;
	else {
		printf("Unrecognized fitting method found \'%s\', should be "
				"\'simplex\' or \'multinest\' or \'nr\'.", 
				fit_method_option->sval[0]);
		return;
	}
	__parsed_ok=true;
}

double FitSubpixCommandLineOptions::siman_option(SIMAN_OPTION option) const
{
	switch(option) {
		case NTRIES: return __siman_ntries;
		case MAX_STEP: return __siman_max_step;
		case BOLTZMAN_K: return __siman_boltzman_k;
		case START_TEMPERATURE: return __siman_start_temperature;
		case COOLING_RATE: return __siman_cooling_rate;
		case MIN_TEMPERATURE: return __siman_min_temperature;
		default: throw Error::InvalidArgument(
						 "FitSubpixCommandLineOptions::siman_option",
						 "Unrecognized simulated annealing options "
						 "encountered");
	}
}


