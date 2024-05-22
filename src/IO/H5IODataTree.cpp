/**\file
 *
 * \brief The definitions of some of the methods of the H5IODataTree class.
 *
 * \ingroup IO
 */

#include "H5IODataTree.h"
#include <iostream>

namespace boost {
  bool operator==(const any& lhs, const any& rhs){
    throw Error::NotImplemented("Comparison of boost any is undefined");
  }
}

namespace IO {

    void H5IODataTree::add_1d_entry(
        double *value,
        unsigned length,
        const std::string &quantity
    )
    {
        if(length > 1) {
            Eigen::VectorXd *entry=new Eigen::VectorXd(length);
            for(unsigned i=0; i<length; ++i)
                (*entry)[i] = value[i];

            put(quantity,
                entry,
                IO::TranslateToAny< Eigen::VectorXd >());
            __doubles_to_destroy.push_back(entry);
        } else if(length == 1) {
            put(quantity,
                value[0],
                IO::TranslateToAny< double >());
        } else
            throw Error::InvalidArgument(
                "add_1d_tree_entry",
                "Attempting to add zero length dataset to I/O tree!"
            );
    }

    void H5IODataTree::initialize_command_line(int argc,
                                               char** argv,
                                               const std::string &executable,
                                               const std::string &version)
    {
        std::ostringstream command_line;
        command_line << "'";
        for(int i = 0; i < argc; ++i)
            command_line << "'" << argv[i] << (i == argc - 1 ? "" : "' ");
        command_line << "'";

        if(executable == "fitpsf" || executable == "fitprf") {
            __tool = (executable == "FitPSF" ? FITPSF : FITPRF);
            __prefix = "psffit.";
            put("bg.tool", executable, translate_string);
            put("bg.cmdline", command_line.str(), translate_string);
            put("bg.version", version, translate_string);
            put("bg.model", "annulus", translate_string);
            put("bg.source_id", "projsrc.srcid.name", translate_string);
            put("bg.source_x", "projsrc.x", translate_string);
            put("bg.source_y", "projsrc.y", translate_string);
        } else if(executable == "subpixphot") {
            __tool = SUBPIXPHOT;
            __prefix = "apphot.";
        } else if(executable == "fitsubpix") {
            __tool = FITSUBPIX;
            __prefix = "fitsubpix.";
        } else assert(false);

        put(__prefix + "tool", executable, translate_string);
        put(__prefix + "cmdline", command_line.str(), translate_string);
        put(__prefix + "version", version, translate_string);
        put(__prefix + "source_id", "projsrc.srcid.name", translate_string);
        put(__prefix + "source_x", "projsrc.x", translate_string);
        put(__prefix + "source_y", "projsrc.y", translate_string);
        put(__prefix + "background", "bg.value", translate_string);
        put(__prefix + "background_err", "bg.error", translate_string);
    }

    void H5IODataTree::process_psffit_option(
        const std::string &key,
        const opt::variable_value &value
    )
    {
        size_t first_split = key.find_first_of('.');
        if(first_split == std::string::npos) {
            std::string destination = key;
            std::replace(destination.begin(), destination.end(), '-', '_');
            put(__prefix + destination, value.value());
        }
        std::string component = key.substr(0, first_split),
                    sub_key = key.substr(first_split + 1);
        std::replace(sub_key.begin(), sub_key.end(), '-', '_');
        if(component == "psf") {
            size_t subkey_split = sub_key.find_first_of('.');
            if(subkey_split == std::string::npos) {
                if(sub_key == "model")
                    put(__prefix + sub_key, __psf_model, translate_string);
                else
                    put(__prefix + sub_key, value.value());
            } else if(sub_key.substr(0, subkey_split) == __psf_model) {
                std::string sub_sub_key = sub_key.substr(subkey_split + 1);
                if(sub_sub_key == "grid")
                    put(__prefix + sub_sub_key,
                        represent_grid(value.as<PSF::Grid>()),
                        translate_string);
                else
                    put(__prefix + sub_sub_key, value.value());
            }
        } else if(component == "io"
                  &&
                  (sub_key == "image"
                   ||
                   sub_key == "source_list"
                   ||
                   sub_key == "initial_guess")
        ) {
            std::string value_str = value.as<std::string>();
            if(value_str != "") {
                put(__prefix + sub_key, value_str, translate_string);
                if(sub_key != "initial_guess")
                    put("bg." + sub_key, value_str, translate_string);
            }
        }
        else if(component == "src")
            put(__prefix + "srcpix_" + sub_key, value.value());
        else if(component == "bg") {
            if(sub_key == "min_pix")
                put(__prefix + "min_bg_pix", value.value());
            else if(sub_key == "annulus")
                put("bg.annulus",
                    value.as<Background::Annulus>(),
                    translate_string);
            else put("bg." + sub_key, value.value());
        }
    }

    void H5IODataTree::process_subpixphot_option(
        const std::string &key,
        const opt::variable_value &value
    )
    {
        if(key == "gain") put(__prefix + "gain", value.value());
        if(key == "io.image") {
            std::string value_str = value.as<std::string>();
            if(value_str != "")
                put("apphot.image", value_str, translate_string);
            return;
        } else if(key == "ap.aperture") {
            Core::RealList aperture_list = value.as<Core::RealList>();
            aperture_list.sort();
            std::valarray<double>
                *apertures = new std::valarray<double>(aperture_list.size());
            unsigned index = 0;
            for(
                Core::RealList::const_iterator ap_i = aperture_list.begin();
                ap_i != aperture_list.end();
                ++ap_i
            )
                (*apertures)[index++] = *ap_i;
            put("apphot.aperture",
                apertures,
                TranslateToAny< std::valarray<double> >());
        } else if(key == "ap.const-error")
            put(__prefix + "const_error", value.value());
        size_t first_split = key.find_first_of('.');
        if(first_split == std::string::npos)
            put(__prefix + key, value.value());
    }

    void H5IODataTree::fill_configuration(const opt::variables_map& options)
    {
        if(__tool==FITPSF || __tool==FITPRF)
            switch(options["psf.model"].as<PSF::ModelType>()) {
                case PSF::SDK : __psf_model="sdk"; break;
                case PSF::BICUBIC : __psf_model="bicubic"; break;
                case PSF::ZERO : __psf_model="zero";
            }
        else __psf_model="";
        for(
                opt::variables_map::const_iterator option_i=options.begin();
                option_i!=options.end();
                ++option_i
        )
            switch(__tool) {
                case FITPSF : case FITPRF :
                    process_psffit_option(option_i->first, option_i->second);
                    break;
                case SUBPIXPHOT :
                    process_subpixphot_option(option_i->first,
                                              option_i->second);
                    break;
                case FITSUBPIX :
                    throw Error::Runtime(
                            "HDF5 I/O for fitsubpix not implemented!"
                    );
            }
    }

    void H5IODataTree::add_1d_string_entry(char **value,
                                           unsigned length,
                                           const std::string &path)
    {
        if(length > 1) {
            std::vector<std::string> *entry = new std::vector<std::string>(
                value,
                value + length
            );
            __strings_to_destroy.push_back(entry);
            put(path,
                entry,
                IO::TranslateToAny< std::vector<std::string> >());
        } else if(length == 1) {
            put(path,
                std::string(value[0]),
                IO::TranslateToAny<std::string>());
        }
        else
            throw Error::InvalidArgument(
                "add_1d_tree_entry",
                "Attempting to add zero length string dataset to I/O tree!"
            );
    }

    void H5IODataTree::add_c_array(const std::string &quantity,
                                   void *value,
                                   const std::string &format,
                                   unsigned length)
    {
        if(format == "str")
            add_1d_string_entry(reinterpret_cast<char**>(value),
                         length,
                         quantity);
        else if(format == "int")
            add_1d_entry(reinterpret_cast<int*>(value),
                         length,
                         quantity,
                         __ints_to_destroy);
        else if(format == "long")
            add_1d_entry(reinterpret_cast<long*>(value),
                         length,
                         quantity,
                         __longs_to_destroy);
        else if(format == "short")
            add_1d_entry(reinterpret_cast<short*>(value),
                         length,
                         quantity,
                         __shorts_to_destroy);
        else if(format == "char")
            add_1d_entry(reinterpret_cast<char*>(value),
                         length,
                         quantity,
                         __chars_to_destroy);
        else if(format == "uint")
            add_1d_entry(reinterpret_cast<unsigned*>(value),
                         length,
                         quantity,
                         __uints_to_destroy);
        else if(format == "ulong")
            add_1d_entry(reinterpret_cast<unsigned long*>(value),
                         length,
                         quantity,
                         __ulongs_to_destroy);
        else if(format == "ushort")
            add_1d_entry(reinterpret_cast<unsigned short*>(value),
                         length,
                         quantity,
                         __ushorts_to_destroy);
        else if(format == "uchar")
            add_1d_entry(reinterpret_cast<unsigned char*>(value),
                         length,
                         quantity,
                         __uchars_to_destroy);
        else if(format == "bool")
            add_1d_entry(reinterpret_cast<bool*>(value),
                         length,
                         quantity,
                         __bools_to_destroy);
        else if(format == "double")
            add_1d_entry(reinterpret_cast<double*>(value),
                         length,
                         quantity);
        else
            throw Error::InvalidArgument(
                "update_result_tree",
                "invalid format: " + std::string(format)
            );

    }



    H5IODataTree::~H5IODataTree()
    {
#ifdef VERBOSE_DEBUG
        std::cerr << "Destroying H5IODataTree at " << this << std::endl;
#endif
        destroy_allocated(__strings_to_destroy);
        destroy_allocated(__ints_to_destroy);
        destroy_allocated(__uints_to_destroy);
        destroy_allocated(__longs_to_destroy);
        destroy_allocated(__ulongs_to_destroy);
        destroy_allocated(__shorts_to_destroy);
        destroy_allocated(__ushorts_to_destroy);
        destroy_allocated(__chars_to_destroy);
        destroy_allocated(__uchars_to_destroy);
        destroy_allocated(__bools_to_destroy);
        destroy_allocated(__doubles_to_destroy);
    }

    std::ostream &operator<<(std::ostream &os, const IOTreeBase &tree)
    {
        static unsigned level = 0;
        const std::string padding = "    ";
        for(
            IOTreeBase::const_iterator node_i = tree.begin();
            node_i != tree.end();
            ++node_i
        ) {
            for(unsigned i = 0; i<level; ++i) os << padding;
            os << "(" << node_i->first << ")"
               << (node_i->second.data().empty() ? "<empty>" : "<filled>")
               << std::endl;
            ++level;
            os << node_i->second;
            --level;
        }
        return os;
    }

}
