/**\file
 * 
 * \brief Defines the function that performs Newton-Raphson fitting for the
 * sub-pixel sensitivity map along with some other functions necessary to
 * accomplish this.
 *
 * \ingroup FitSubpix
 */

#include "NR_fitting.h"

///\brief Constructs a vector of initial sensitivities
///
///Each call returns a different initial map. The first one or several per
///sub-pixel map resolution are pre-defined, after that, random maps are
///returned.
Eigen::VectorXd initial_x(unsigned num_subpix)
{
/*	Eigen::VectorXd s(num_subpix);
	s << 1, 1, 0, 1,
	     1, 0, 0, 1,
		 0, 1, 0, 1,
		 0, 0, 0, 1;
	return s_to_x(s/s.sum());*/

/*	Eigen::VectorXd s(num_subpix);
	s << 0.157131,	0.988891,	0.7807,		0.709905,
	     2.04698,	3.8153,		0.625146,	2.12874,
		 0.616719,	0.440388,	0.246857,	0.909662,
		 1.37552,	0.430847,	0.279564,	0.447651;
	return s_to_x(s/s.sum());*/

/*	Eigen::VectorXd s(num_subpix);
	s << 1.5508,	1.5468,		3.47493,	0.201079,
	     2.9255,	0.124357,	0.572004,	0.312847,
		 0.00852634,2.69834,	0.15837,	0.023461,
		 0.424276,	0.208796,	0.503905,	1.26601;
	return s_to_x(s/s.sum());*/

	static int call_count=0;

	Eigen::VectorXd s(num_subpix);
	if(call_count==0) 
		s=Eigen::VectorXd::Ones(num_subpix);
	else if(call_count==1 && num_subpix==16)
		s << 2.69542e-01, 2.69542e-01, 1.13882e+00, 1.13882e+00,
		     2.69542e-01, 2.69542e-01, 1.13882e+00, 1.13882e+00,
			 1.75950e+00, 1.75950e+00, 8.32139e-01, 8.32139e-01,
			 1.75950e+00, 1.75950e+00, 8.32139e-01, 8.32139e-01;
	else if(call_count==1 && num_subpix==9)
		s << 2.69542e-01, 0.5*(2.69542e-01+1.13882e+00), 1.13882e+00,
		     0.5*(2.69542e-01+1.75950e+00),	0.25*(2.69542e-01+1.13882e+00+1.75950e+00+8.32139e-01), 0.5*(1.13882e+00+8.32139e-01),
			 1.75950e+00, 0.5*(1.75950e+00+8.32139e-01), 8.32139e-01;
	else if(call_count==1 && num_subpix==25) {
		s << 2.59341e+00, 2.77012e-01, 5.52777e-01, 1.32756e-01, 3.51339e+00,
		     6.01402e-01, 0.00000e+00, 3.87031e-02, 0.00000e+00, 0.00000e+00,
			 5.09771e-02, 5.26300e-01, 0.00000e+00, 0.00000e+00, 4.87233e-01,
			 4.10883e-01, 4.76827e+00, 2.13416e+00, 1.89672e+00, 4.04353e-01,
			 0.00000e+00, 4.15542e+00, 1.08776e+00, 6.74921e-01, 6.93551e-01;
		std::cout << "Checking siman subpix" << std::endl;
	} else if(call_count==1 && num_subpix==36)
		s << 2.50225e-01, 5.31725e-02, 2.24446e+00, 3.02764e+00, 2.79883e-01,
			 												0.00000e+00,
			 1.86050e-01, 1.46305e-01, 1.56628e+00, 9.35313e-02, 2.32570e+00,
			 												7.27272e-02,
			 3.89608e-01, 0.00000e+00, 2.49588e-01, 8.89488e-01, 9.75034e-01,
			 												9.66325e-01,
			 1.26314e+00, 9.25585e-01, 8.99482e-01, 0.00000e+00, 9.57012e-01,
			 												3.42674e-01,
			 3.33602e-01, 3.17664e+00, 0.00000e+00, 0.00000e+00, 1.87037e+00,
			 												2.29400e-01,
			 6.52785e+00, 2.47943e-01, 8.14227e-01, 2.74401e-01, 2.41400e+00,
		  													2.00766e+00;
	else if(call_count==2 && num_subpix==36)
		s << 0.00000e+00, 2.21549e+00, 1.74378e+00, 0.00000e+00, 2.62744e-01,
			 												6.01035e-01,
			 1.44025e-01, 3.08675e-01, 1.37542e+00, 1.35779e+00, 1.93066e+00,
			 												0.00000e+00,
			 0.00000e+00, 4.67434e-02, 1.82137e-01, 1.28451e+00, 6.72736e-01,
			 												0.00000e+00,
		 	 2.90404e-01, 2.91433e+00, 8.83459e-01, 1.58258e+00, 3.16820e-01,
			 												0.00000e+00,
			 3.16544e-01, 1.60408e+00, 5.65125e+00, 9.79755e-01, 1.62900e+00,
			 												1.21893e+00,
			 2.59919e+00, 5.71511e-01, 7.64609e-01, 3.99847e-01, 1.58110e+00,
		  													5.70836e-01;
	else if(call_count==2 && num_subpix==16)
		s << 2.37406e-01,
		     0.25*2.37406e-01+0.75*2.63276e+00,
			 0.75*2.63276e+00+0.25*9.57954e-03,
			 9.57954e-03,

		     0.25*2.37406e-01+0.75*0,
		     0.25*(0.25*2.37406e-01+0.75*2.63276e+00)+0.75*(0.25*0+0.75*1.43128e-01),
			 0.25*(0.75*2.63276e+00+0.25*9.57954e-03)+0.75*(0.75*1.43128e-01+0.25*1.07927e+00),
			 0.25*9.57954e-03+0.75*1.07927e+00,

			 0.25*0+0.75*3.13010e+00,
		     0.25*(0.25*0+0.75*1.43128e-01)+0.75*(0.25*3.13010e+00+0.75*1.96497e-01),
			 0.25*(0.75*1.43128e-01+0.25*1.07927e+00)+0.75*(0.25*1.96497e-01+0.75*1.57126e+00),
			 0.25*1.07927e+00+0.75*1.57126e+00,
			 
			 3.13010e+00,
		     0.25*3.13010e+00+0.75*1.96497e-01,
			 0.75*1.96497e-01+0.25*1.57126e+00,
			 1.57126e+00;
	else if(call_count==3 && num_subpix==16)
		s << 1, 1,  1,  1,
		     1, 10, 10, 1,
			 1, 10, 10, 1,
			 1, 1,  1,  1;
	else 
		for(unsigned i=0; i<num_subpix; i++) 
			s(i)=static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	call_count++;
	return s_to_x(num_subpix*s/s.sum());

/*	int N=num_subpix-1;
	Eigen::VectorXd result(N);
	const double s=1.0/num_subpix;
	for(int i=0; i<N; i++) 
		result(i)=1.0-std::pow(1.0-s/(1.0-s*i), N-i);
	return result;*/
}

///\brief Calcualtes by how much dx should be scaled so that at the updated
///point the variance is less than at the present point.
///
///Also updates var_val with the variance at the new point.
double dx_scaling(const Eigen::VectorXd &x, const Eigen::VectorXd &dx, 
		Variance &var, double &var_val, int direction=1)
{
	const double min_dx=1e-5;
	double gamma=1.0, diff=1.0, 
		   min_gamma=min_dx/dx.norm();
	unsigned N=dx.size();
	double var_val_new = 0;
	while(diff>0 && gamma>min_gamma) {
		gamma/=2.0;
		for(unsigned i=0; i<N; i++) {
			double dxi=direction*dx(i);
			if(dxi>0) gamma=std::min(gamma, 0.5*(1.0-x(i))/dxi);
			else if(dxi<0) gamma=std::min(gamma, -x(i)/dxi);
		}
		var_val_new=var.variance(x_to_s(x+direction*gamma*dx));
		diff=var_val_new-var_val;
	}
	if(gamma<=min_gamma) {
		if(direction==-1) return 0;
		else return dx_scaling(x, dx, var, var_val, -1);
	} else {
		var_val=var_val_new;
		return direction*gamma;
	}
}

///\brief Updates x according to a Newton-Raphson step and returns the norm
///of the change.
double update_x(Eigen::VectorXd &x, Variance &var, double &var_val,
		double var_ref)
{
	unsigned num_subpix=x.size()+1;
	Eigen::VectorXd var_grad, s=x_to_s(x);
	Eigen::MatrixXd var_d2, jacobian=calculate_jacobian(x);
	var.calculate(s, &var_val, &var_grad, &var_d2);
	Eigen::VectorXd RHS=-jacobian*var_grad;
	Eigen::MatrixXd d2=jacobian*var_d2*jacobian.transpose();
	std::valarray<Eigen::MatrixXd> d2s_dx2=calculate_d2s(x); 
	for(unsigned i=0; i<num_subpix; i++) d2+=var_grad(i)*d2s_dx2[i];
	Eigen::FullPivHouseholderQR<Eigen::MatrixXd> decomposition(d2);
	Eigen::VectorXd dx=decomposition.solve(RHS);
	double gamma=dx_scaling(x, dx, var, var_val);
	std::cout << "variance=" << var_val << ", " << var_ref/var_val
		<< " times better than reference." << std::endl;
	dx*=gamma;
	x+=dx;//.cwiseMax(-x);
	return dx.norm();
}

void fit_using_NR(const FitSubpixCommandLineOptions &options)
{
	const double tolerance=1e-8;
	Variance var(options.frames(), options.sources(), options.aperture(),
			options.input_columns(), options.background_annulus(), 
			options.x_split(), options.y_split());
	SubPixelMap best_map("Best map");
	double best_variance=Inf, uniform_variance;
	int best_run = -1;
	srand(time(NULL));
	for(int run=0; run<10; run++) {
		Eigen::VectorXd x=initial_x(options.x_split()*options.y_split());
		Eigen::VectorXd s=x_to_s(x);
		if(run==0) {
			var.calculate(s, &uniform_variance);
			std::cout << "Variance for uniform subpixel sensitivity: "
				<< uniform_variance << " used as reference." << std::endl;
		}

		SubPixelMap map(options.x_split(), options.y_split(), "solution");
		for(unsigned y=0; y<options.y_split(); y++) 
			for(unsigned x=0; x<options.x_split(); x++)
				map(x,y)=s(x+y*options.x_split());
		std::cout << "Run " << run << " initial sub-pixel map:" << std::endl
			<< map << std::endl;
		double dx_norm=tolerance+1, final_variance;
		while(dx_norm>tolerance) dx_norm=update_x(x, var, final_variance,
				uniform_variance);
		s=x_to_s(x);
		for(unsigned y=0; y<options.y_split(); y++) 
			for(unsigned x=0; x<options.x_split(); x++)
				map(x,y)=s(x+y*options.x_split());
		std::cout << "Run " << run << " converged to variance=" 
			<< final_variance << ", for a subpixel map:" 
			<< std::endl << map << std::endl;
		if(final_variance<best_variance) {
			best_variance=final_variance;
			best_map=map;
			best_run=run;
		}
	}
	std::cout << "The best subpixel map was found on run" << best_run 
		<<": " << std::endl << best_map
		<< " and it leads to a variance of: " << best_variance << std::endl;
}
