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

#include "../Core/SharedLibraryExportMacros.h"
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
	class LIB_PUBLIC General : public std::exception {
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
    class LIB_PUBLIC NotImplemented : public General {
    public:
        NotImplemented(const std::string &error_message="")
            : General(error_message)
            {}
		virtual const char *what() const throw()
        {return "Not implemented error";}
    };

	///\brief The parent of all fits library related errors.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC Fits : public General {
	public:
		Fits(const std::string &error_message="") : General(error_message) {}
		virtual const char *what() const throw() {return "Fits file error";}
	};

	///\brief Errors related an the actual fits image.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC FitsImage : public Fits {
	public:
		FitsImage(const std::string &error_message="") : Fits(error_message) {}
		virtual const char *what() const throw () {return "Fits image error";}
	};

	///\brief An attempt was made to access a pixel outside the image area.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC ImageOutside : public FitsImage {
	public:
		ImageOutside(unsigned long x,
                     unsigned long y,
                     unsigned long xres,
                     unsigned long yres)
		{
			std::ostringstream msg;
			msg << "Attempting to access outside the image area ("
				<< xres << " by " << yres << "): " << x << ", " << y;
			set_message(msg.str());
		}
		virtual const char *what() const throw() {return "Outside fits image error";}
	};

	///\brief The parent of all run-time errors.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC Runtime : public General {
	public:
		Runtime(const std::string &error_message="") : General(error_message) {}
		virtual const char *what() const throw() {return "Runtime error";}
	};

	///\brief The type of something was not what was expected.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC Type : public Runtime {
	public:
		Type(const std::string &error_message="") : Runtime(error_message) {}
		virtual const char *what() const throw()
		{return "Unexpected type error";}
	};

	///\brief A function (or method) received an argument with an invalid value.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC InvalidArgument : public Runtime {
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
	class LIB_PUBLIC CommandLine : public Runtime {
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
    class LIB_PUBLIC ParsingError : public Runtime {
    public:
        ParsingError(const std::string &error_message = "") :
            Runtime(error_message) {}
		virtual const char *what() const throw() {return "Bad expression";}
    };

	///\brief Input/Output error.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC IO : public Runtime {
	public:
		IO(const std::string &error_message="") :
				Runtime(error_message) {}
		virtual const char *what() const throw()
		{return "Failed I/O operation";}
	};

	///\brief %Error in a <a href="http://www.gnu.org/software/gsl/">
	///GNU Scientific Library function.</a>
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC GSLError : public Runtime {
	public:
		GSLError(const std::string &error_message="") :
			Runtime(error_message) {}
		virtual const char *what() const throw() {return "GSL Error";}
	};

	///%Error while fitting.
	///\ingroup SubPixPhot
	///\ingroup FitSubpix
	///\ingroup FitPSF
	class LIB_PUBLIC Fitting : public Runtime {
	public:
		Fitting(const std::string &error_message="") :
			Runtime(error_message) {}
		virtual const char *what() const throw() {return "Fitting failed";}
	};
};

#endif
