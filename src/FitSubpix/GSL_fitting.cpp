/**\file
 *
 * \brief Implementation of the simplex and simulated annealing fitting for
 * the sub-pixel structure.
 *
 * \ingroup FitSubpix
 */

#include "GSL_fitting.h"

void GSLSubPixelMap::set_sensitivities(const gsl_vector *cube)
{
	unsigned long x_res=x_resolution(), y_res=y_resolution(), 
				  nsubpix=x_res*y_res, N=nsubpix-1;
	if(cube->size<nsubpix-1)
		throw Error::InvalidArgument("GSLSubPixelMap constructor", 
				"Invalid size of cube vector, should be at least "
				"the number of subpixels minus one");
	double range=1.0;
	for(unsigned i=0; i<nsubpix-1; i++) {
		unsigned long x=i%x_res, y=i/x_res;
		double cube_val=gsl_vector_get(cube, i); 
		if(cube_val<0.0 || cube_val>1.0) {
			std::ostringstream arg_msg;
			arg_msg.precision(16);
			arg_msg << "cube[" << i << "]=" << cube_val;
			throw Error::InvalidArgument(
				"GSLSubPixelMap::set_sensitivities", arg_msg.str());
		}
		double sensitivity=range*(1.0-std::pow(1.0-cube_val, 1.0/(N-i)));
		range-=sensitivity;
		(*this)(x,y)=sensitivity;
	}
	(*this)(x_res-1, y_res-1)=range;
}

///The function to be minimized by GSL
double GSLchi2_func(const gsl_vector *fit_var, void *params)
{
	void **param_array=static_cast<void **>(params);
	ChiSquared *chi2=static_cast<ChiSquared *>(param_array[0]);
	unsigned xres=*(static_cast<unsigned *>(param_array[1])), 
		 yres=*(static_cast<unsigned *>(param_array[2]));
	static GSLSubPixelMap subpix_map(xres, yres);
	try {subpix_map.set_sensitivities(fit_var);}
	catch (Error::InvalidArgument) {
		gsl_vector *closest_inside=gsl_vector_alloc(fit_var->size);
		double distance=0.0;
		double eps=std::numeric_limits<double>::epsilon();
		std::cout.precision(16);
		for(size_t i=0; i<fit_var->size; i++) {
			double value=gsl_vector_get(fit_var, i);
			if(value>1.0) {
				gsl_vector_set(closest_inside, i, 1.0);
				distance+=(value-1.0)/eps;
			} else if(value<0.0) {
				gsl_vector_set(closest_inside, i, 0.0);
				distance+=-value/eps;
			} else gsl_vector_set(closest_inside, i, value);
		}
		subpix_map.set_sensitivities(closest_inside);
		return (*static_cast<ChiSquared *>(params))(subpix_map)*
				std::exp(distance);
	}
	if(fit_var->size==xres*yres+4)
		chi2->set_SDK_amp(1.0/gsl_vector_get(fit_var, xres*yres),
						1.0-2.0*gsl_vector_get(fit_var, xres*yres+1),
						1.0-2.0*gsl_vector_get(fit_var, xres*yres+2),
						gsl_vector_get(fit_var, xres*yres+3));
	else if(fit_var->size==xres*yres+1)
		chi2->set_SDK_amp(gsl_vector_get(fit_var, xres*yres),
						1.0/gsl_vector_get(fit_var, xres*yres), 0.0, 0.0);
	else if(fit_var->size!=xres*yres) {
		std::ostringstream msg;
		msg << "Unexpected number of fit parameters ("
			<< fit_var->size << "). Should be either "
			<< xres*yres-1 << " or " << xres*yres+3;
		throw Error::InvalidArgument("GSLchi2_func", msg.str());
	}
	double result=(*chi2)(subpix_map);
	std::cout << "Chi2(";
	for(size_t i=0; i<fit_var->size-1; i++)
		std::cout << gsl_vector_get(fit_var, i) << ", ";
	std::cout << gsl_vector_get(fit_var, fit_var->size-1) 
			<< ")=" << result << std::endl;
	return result;
}

void fit_using_GSL_simplex(const FitSubPixConfig &options)
{
	unsigned dimensions=options["subpix.x-split"].as<unsigned>
                        *
                        options["subpix.y_split"].as<unsigned>-1;

	ChiSquared chi2(options["io.image-series"].as<StringList>(), 
                    options["io.psffit-series"].as<StringList>(), 
                    options["io.subpix.aperture"].as<double>(),
					options["psf.sdk.max_exp_coef"].as<double>());

	gsl_multimin_function minimizer_func;

	gsl_vector *initial_sensitivities=gsl_vector_alloc(dimensions);
	gsl_vector_set_all(initial_sensitivities, 0.0);

	minimizer_func.n = (dimensions);
	minimizer_func.f = GSLchi2_func;
	unsigned xres=options.x_split(), yres=options.y_split();
	void *func_params[]={&chi2, &xres, &yres};
	minimizer_func.params=func_params;

	gsl_vector *step_sizes = gsl_vector_alloc(dimensions);
	gsl_vector_set_all(step_sizes, 1.0);

	gsl_multimin_fminimizer *minimizer = gsl_multimin_fminimizer_alloc(
					gsl_multimin_fminimizer_nmsimplex2, dimensions); 
	gsl_multimin_fminimizer_set(minimizer, &minimizer_func, 
					initial_sensitivities, step_sizes);
	int status=GSL_CONTINUE;
	for(unsigned iter=0; status==GSL_CONTINUE; iter++) {
		status = gsl_multimin_fminimizer_iterate(minimizer);
		if (status) break;

		double size = gsl_multimin_fminimizer_size(minimizer);
		status = gsl_multimin_test_size (size, 1e-4);
		if(status == GSL_SUCCESS) 
			std::cout << "converged to minimum at" << std::endl;
		std::cout << std::setw(4) << iter << " " << minimizer << std::endl;
		status=GSL_CONTINUE;
	}
	
	gsl_vector_free(initial_sensitivities);
	gsl_vector_free(step_sizes);
	gsl_multimin_fminimizer_free (minimizer);
}

///The energy to pass to the simulated annealing minimizer of GSL
double simulated_annealing_energy(void *configuration)
{
	void **var_x=static_cast<void **>(configuration);
	Variance *var=static_cast<Variance *>(var_x[0]);
	Eigen::VectorXd *x=static_cast<Eigen::VectorXd *>(var_x[1]);
	double result;
	var->calculate(x_to_s(*x), &result);
	return result;
}

///\brief Distance between two unit-cube representations of sub-pixel
///sensitivity maps (usual L2 norm).
double simulated_annealing_distance(void *config1, void *config2)
{
	Eigen::VectorXd
		*x1=static_cast<Eigen::VectorXd *>(static_cast<void **>(config1)[1]),
		*x2=static_cast<Eigen::VectorXd *>(static_cast<void **>(config2)[1]);
	return (*x1-*x2).norm();
}

///The stepping function for the simulated annealing minimizer
void simulated_annealing_step(const gsl_rng * rand_number_generator,
		void *configuration, double step_size)
{
	Eigen::VectorXd *x=static_cast<Eigen::VectorXd *>(
			static_cast<void **>(configuration)[1]);
	Eigen::VectorXd change(x->size());
	for(int i=0; i<x->size(); i++) {
		change[i] = 2.0*gsl_rng_uniform(rand_number_generator)-1.0;
	}
	change*=step_size/change.norm();
	*x+=change;
	for(int i=0; i<x->size(); i++)
		if((*x)[i]<0) (*x)[i]=0;
		else if((*x)[i]>1) (*x)[i]=1;
}

///\brief Copies GSL simulated annealing configuration source to dest.
///
///Copies the subpixel structure, but leaves the variance variable alone.
void simulated_annealing_copy_config(void *source, void *dest)
{
	*static_cast<Eigen::VectorXd *>(static_cast<void **>(dest)[1])=
		*static_cast<Eigen::VectorXd *>(static_cast<void **>(source)[1]);
}

///Creates a copy of a GSL simulated annealing configuration.
void *simulated_annealing_copy_construct(void *configuration)
{
	void **var_x=static_cast<void **>(configuration);
	void **result=new void*[2];
	result[0]=var_x[0];
	result[1]=static_cast<void *>(new Eigen::VectorXd(
				*static_cast<Eigen::VectorXd*>(var_x[1])));
	return static_cast<void*>(result);
}

///\brief Destroys a GSL simulated annealing configuration.
///
///Deletes the subpixel part, leaving the variance variable alone.
void simulated_annealing_destroy(void *configuration)
{
	delete static_cast<Eigen::VectorXd *>(
			static_cast<void **>(configuration)[1]);
}

///Output function for the simulated annealing configuration.
void simulated_annealing_output(void *configuration)
{
	std::cout << (*static_cast<Eigen::VectorXd *>(
			static_cast<void **>(configuration)[1])).transpose();
	std::cout.flush();
}

void fit_using_GSL_simulated_annealing(const FitSubPixConfig &options)
{
    Variance var(
        options.frames(),
        options.sources(),
        options.aperture(),
        options.input_columns(),
        options.background_annulus(),
        options.x_split(),
        options.y_split()
    );
	unsigned num_subpix=options.x_split()*options.y_split();
	void *initial_config[2];
	initial_config[0]=static_cast<void*>(&var);
	initial_config[1]=static_cast<void*>(
			new Eigen::VectorXd(s_to_x(Eigen::VectorXd::Ones(num_subpix))));

	gsl_siman_params_t simulated_annealing_parameters 
		= {static_cast<int>(options.siman_option(NTRIES)),
			1,	//Decrease the temperature every iterations (simply use 
			    //slower cooling rate if required).
			options.siman_option(MAX_STEP),
			options.siman_option(BOLTZMAN_K),
			options.siman_option(START_TEMPERATURE),
			options.siman_option(COOLING_RATE),
			options.siman_option(MIN_TEMPERATURE)};

	const gsl_rng_type * RNGType;
	gsl_rng * random_number_generator;
     
	gsl_rng_env_setup();
     
	RNGType = gsl_rng_default;
	random_number_generator = gsl_rng_alloc(RNGType);
     
	gsl_siman_solve(random_number_generator,
			static_cast<void*>(initial_config),
			simulated_annealing_energy, simulated_annealing_step,
			simulated_annealing_distance, simulated_annealing_output,
			simulated_annealing_copy_config,
			simulated_annealing_copy_construct, simulated_annealing_destroy,
			0, simulated_annealing_parameters);
     
	gsl_rng_free (random_number_generator);
}

std::ostream &operator<<(std::ostream &os, 
				gsl_multimin_fminimizer *minimizer)
{
	os.precision(4);
	static GSLSubPixelMap subpix_map(4,4);
	try {
		subpix_map.set_sensitivities(minimizer->x);
		os << subpix_map;
	} catch (Error::InvalidArgument) {}
/*	for(size_t i=0; i<minimizer->x->size; i++) {
		if(minimizer) os << std::setw(7) << gsl_vector_get(minimizer->x, i);
		else os << std::setw(3) << "s[" << std::setw(3) << i << "]";
		os << " ";
	}*/
	return os;
}

