/**\file
 *
 * \brief The hierarchy of exceptions for this project.
 *
 * All exceptions have the usual what() member function as well as an
 * additional get_message() which provides more information on what went
 * wrong.
 *
 * \ingroup SubPixPhot
 * \ingroup FitSubpix
 * \ingroup FitPSF
 */

#ifndef __ERROR_H
#define __ERROR_H

#include "fitsio.h"
#include <sstream>
#include <iostream>

/**\brief The namespace containing the exception hierarchy.
 *
 * \ingroup SubPixPhot
 * \ingroup FitSubpix
 * \ingroup FitPSF
 */
namespace Error {
	///\brief The parent of all errors raised anywhere in the code.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class General : public std::exception {
	private:
		std::string message;
	public:
		General(const std::string &error_message="") : 
			message(error_message) 
            {
#ifndef NDEBUG
                std::cerr << what() << ": " << get_message() << std::endl;
#endif
            }
		virtual void set_message(const std::string &error_message)
		{message=error_message;}
		virtual const char *what() const throw() {return "General error";}
		virtual const std::string &get_message() {return message;}
		virtual ~General() throw() {}
	};

    ///\brief An error indicating a feature has not been implemented yet.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
    class NotImplemented : public General {
    public:
        NotImplemented(const std::string &error_message="")
            : General(error_message)
            {}
		virtual const char *what() const throw() {return "Fits file error";}
    };

	///\brief The parent of all fits library related errors.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class Fits : public General {
	public:
		Fits(const std::string &error_message="") : General(error_message) {}
		virtual const char *what() const throw() {return "Fits file error";}
	};

	///\brief Errors related an the actual fits image.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class FitsImage : public Fits {
	public:
		FitsImage(const std::string &error_message="") : Fits(error_message) {}
		virtual const char *what() const throw () {return "Fits image error";}
	};

	///\brief An attempt was made to access a pixel outside a fits image area.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class FitsImageOutside : public FitsImage {
	public:
		FitsImageOutside(unsigned long x, unsigned long y, unsigned long xres,
				unsigned long yres, 
				const std::string &filename="unknown fits file", 
				int hdu_number=0) 
		{
			std::ostringstream msg;
			msg << "Attempting to access outside the image area ("
				<< xres << " by " << yres << "of HDU ";
			if(hdu_number) msg << hdu_number;
			else msg << "unknown";
			msg << " of '" << filename << "': " << x << ", " << y;
			set_message(msg.str());
		}
		virtual const char *what() const throw() {return "Outside fits image error";}
	};
	
	///\brief The parent of all run-time errors.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class Runtime : public General {
	public:
		Runtime(const std::string &error_message="") : General(error_message) {}
		virtual const char *what() const throw() {return "Runtime error";}
	};

	///\brief The type of something was not what was expected.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class Type : public Runtime {
	public:
		Type(const std::string &error_message="") : Runtime(error_message) {}
		virtual const char *what() const throw()
		{return "Unexpected type error";}
	};

	///\brief A function (or method) received an argument with an invalid value.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class InvalidArgument : public Runtime {
	public:
		InvalidArgument(
				///Name of the function that received the invalid argument.
				const std::string &func_name,
				
				///A message giving details about which argument and why.
				const std::string &arg_msg)
		{
			std::ostringstream msg;
			msg << arg_msg << " in " << func_name;
			set_message(msg.str());
		}

		InvalidArgument(
				///Name of the function that received the invalid argument.
				const std::string &func_name,
				
				///The number of the offending argument.
				int arg_num)
		{
			std::ostringstream msg;
			msg << "argument " << arg_num << " to " << func_name;
			set_message(msg.str());
		}
		virtual const char *what() const throw()
		{return "Invalid function argument";}
	};

	///\brief %Error while parsing the command line.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class CommandLine : public Runtime {
	public:
		CommandLine(const std::string &error_message="") : 
				Runtime(error_message) {}
		virtual const char *what() const throw()
        {return "Bad command line or config";}
	};

    ///\brief %Error while parsing the configuration.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
    class ParsingError : public Runtime {
    public:
        ParsingError(const std::string &error_message = "") :
            Runtime(error_message) {}
		virtual const char *what() const throw() {return "Bad expression";}
    };

	///\brief Input/Output error.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class IO : public Runtime {
	public:
		IO(const std::string &error_message="") : 
				Runtime(error_message) {}
		virtual const char *what() const throw()
		{return "Failed I/O operation";}
	};

	///\brief Input/Output error with an HDF5 file.
	class HDF5 : public IO {
	private:
		std::string __path, __message;
	public:
		HDF5(	
				///The path (including the filename) where the error
				///occurred.
				const std::string &path="",

				///Message about what went wrong.
				const std::string &error_message="") : 
			IO(error_message + " at " + path), 
			__path(path),
			__message(error_message) 
            {}

		virtual const char *what() const throw()
		{return "Failed HDF5 file I/O operation:";}

		const std::string &get_path() const {return __path;}
		void set_path(const std::string &path)
		{
			__path=path;
			set_message(__message + " at " + path);
		}

		virtual ~HDF5() throw() {}
	};

	///Signal that an expected component does not exist.
	class HDF5NotFound : public HDF5 {
	public:
		HDF5NotFound(
				///The path (including the filename) where the error
				///occurred.
				const std::string &path,

				///The key this component is identified by.
				const std::string &key,

				///The component missing (e.g. "group" or "dataset" or ...)
				const std::string &component_type) :
			HDF5(path, component_type + "with key '" + key + "' not found!")
			{}

		virtual const char *what() const throw()
		{return "Component missing from HDF5 file.";}
	};

	///\brief %Error in a <a href="http://www.gnu.org/software/gsl/">
	///GNU Scientific Library function.</a>
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class GSLError : public Runtime {
	public:
		GSLError(const std::string &error_message="") :
			Runtime(error_message) {}
		virtual const char *what() const throw() {return "GSL Error";}
	};

	///\brief %Error related to the
	///<a href="http://heasarc.gsfc.nasa.gov/fitsio/">
	///CFITSIO.</a> library.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class CFITSIO : public Runtime {
	public:
		CFITSIO(int cfitsio_error_code,
				const std::string &error_message="")
		{
			std::ostringstream msg;
			char cfitsio_msg[81];
			fits_get_errstatus(cfitsio_error_code, cfitsio_msg);
			msg << error_message << "CFITSIO Error code "
				<< cfitsio_error_code << ": " << cfitsio_msg;
			while(fits_read_errmsg(cfitsio_msg))
				msg << ": " << cfitsio_msg;
			set_message(msg.str());
		}
		virtual const char *what() const throw() {return "CFITSIO Error";}
	};

	///%Error while fitting.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class Fitting : public Runtime {
	public:
		Fitting(const std::string &error_message="") :
			Runtime(error_message) {}
		virtual const char *what() const throw() {return "Fitting failed";}
	};
};

#endif
