/**\file
 *
 * \brief The definitions of some of the methods of the SubPixHDF5File class.
 *
 * \ingroup IO
 */

#include "SubPixHDF5File.h"

namespace IO {

    boost::property_tree::ptree	SubPixHDF5File::__structure,
                                SubPixHDF5File::__dataset_paths,
                                SubPixHDF5File::__attribute_paths,
                                SubPixHDF5File::__link_paths;

    const H5::StrType SubPixHDF5File::VARLEN_STR_TYPE(H5::PredType::C_S1,
                                                      H5T_VARIABLE);

    const hsize_t SubPixHDF5File::__chunk_size=100;

    template<typename SCALAR_TYPE>
    std::ostream& operator<<(std::ostream &os,
                             const OutputArray<SCALAR_TYPE> &array)
    {
        for(unsigned i=0; i<array.size(); ++i)
            os << array[i] << std::endl;
        return os;
    }

    void simple_read_from_h5(const H5::DataSet &dataset,
                             const H5::DataType &memory_type,
                             void *destination)
    {
        dataset.read(destination, memory_type);
    }

    void simple_read_from_h5(const H5::Attribute &attribute,
                             const H5::DataType &memory_type,
                             void *destination)
    {
        attribute.read(memory_type, destination);
    }

    void output_tree(const boost::property_tree::ptree &pt,
            const std::string &indent="")
    {
        using boost::property_tree::ptree;
        for(ptree::const_iterator node=pt.begin(); node!=pt.end(); ++node) {
            std::cout << indent << node->first << ": "
                      << node->second.data() << std::endl;
            output_tree(node->second, indent+"    ");
        }
    }

    SubPixHDF5File::__H5COMPONENT SubPixHDF5File::get_node_type(
            const boost::property_tree::ptree &node)
    {
        std::string type=node.get<std::string>("<xmlattr>.type");
        if( 	type=="group")		return GROUP;
        else if(type=="dataset")	return DATASET;
        else if(type=="link")		return LINK;
        else if(type=="attribute")	return ATTRIBUTE;
        else
            throw Error::Type("Invalid type: " + type + " encountered in "
                              + "SubPixHDF5File::node_type()");
    }

    const H5::DataType &SubPixHDF5File::dataset_type(const std::string &dataset_key)
    {
        if(dataset_key == "projsrc.srcid.field")
            return H5::PredType::STD_U16LE;
        else if(dataset_key == "projsrc.srcid.source") 
            return H5::PredType::STD_U32LE;
        else if(dataset_key == "projsrc.srcid.name")
            return VARLEN_STR_TYPE;
        else if(
                dataset_key	==	"projsrc.x"
                || 	dataset_key	==	"projsrc.y"
                ||	dataset_key	==	"bg.value"
                ||	dataset_key	==	"bg.error"
                ||	dataset_key	==	"psffit.subpixmap"
                ||	dataset_key	==	"psffit.psfmap"
                ||	dataset_key	==	"psffit.variables"
                ||	dataset_key	==	"psffit.mag"
                ||	dataset_key	==	"psffit.flux"
                ||	dataset_key	==	"psffit.amplitude"
                ||	dataset_key	==	"psffit.mag_err"
                ||	dataset_key	==	"psffit.flux_err"
                ||	dataset_key	==	"psffit.amplitude_err"
                ||	dataset_key	==	"psffit.mask_mag"
                ||	dataset_key	==	"psffit.mask_flux"
                ||	dataset_key	==	"psffit.mask_amplitude"
                ||	dataset_key	==	"psffit.mask_mag_err"
                ||	dataset_key	==	"psffit.mask_flux_err"
                ||	dataset_key	==	"psffit.mask_amplitude_err"
                ||	dataset_key	==	"psffit.chi2"
                ||	dataset_key	==	"psffit.sigtonoise"
                ||	dataset_key	==	"psffit.s"
                ||	dataset_key	==	"psffit.d"
                ||	dataset_key	==	"psffit.k"
                ||	dataset_key	==	"apphot.subpixmap"
                ||	dataset_key	==	"apphot.mag"
                ||	dataset_key	==	"apphot.flux"
                ||	dataset_key	==	"apphot.mag_err"
                ||	dataset_key	==	"apphot.flux_err"
        )
            return H5::PredType::IEEE_F32LE;
        else if(
                dataset_key	==	"psffit.npix"
                ||	dataset_key	==	"bg.npix"
        )
            return H5::PredType::STD_U32LE;
        else if(
                dataset_key	==	"psffit.quality"
                ||	dataset_key	==	"apphot.quality"
        )
            return H5::PredType::STD_U8LE;
        else throw Error::InvalidArgument(
                "SubPixHDF5File::dataset_type()",
                "Urecognized dataset key '" + dataset_key +"'"
        );
    }

    void SubPixHDF5File::fill_paths(
            boost::property_tree::ptree& sub_structure,
            const std::string& path,
            const __H5COMPONENT parent_type)
    {
        using boost::property_tree::ptree;
        for(ptree::iterator node=sub_structure.begin();
                node!=sub_structure.end(); ++node) {
            if(node->first!="<xmlattr>") {
                std::string element_name=node->second.get(
                        "<xmlattr>.element_name",
                        node->first
                );
                std::string node_path=(path=="" ? "" : path+"/")
                                      + 
                                      element_name;

                std::string type=node->second.get<std::string>(
                    "<xmlattr>.type"
                );
                if(	type=="dataset") {
                    std::string dataset_id=node->second.get<std::string>(
                            "<xmlattr>.value"
                    );
                    __dataset_paths.put(dataset_id, node_path);
                    if(dataset_id=="projsrc.srcid") {
                        dataset_id="projsrc.srcid.name";
                        node->second.put("<xmlattr>.value", dataset_id);
                        node->second.put("field.<xmlattr>.type", "dataset");
                        node->second.put("field.<xmlattr>.value",
                                         "projsrc.srcid.field");
                        node->second.put(
                                "field.<xmlattr>.compression",
                                node->second.get<std::string>(
                                    "<xmlattr>.compression",
                                    ""
                                )
                        );
                        node->second.put("source.<xmlattr>.type", "dataset");
                        node->second.put("source.<xmlattr>.value",
                                         "projsrc.srcid.source");
                        node->second.put(
                                "source.<xmlattr>.compression",
                                node->second.get<std::string>(
                                    "<xmlattr>.compression",
                                    ""
                                )
                        );
                    } 
                } else if(type=="attribute") {
                    __attribute_paths.put(
                        node->second.get<std::string>("<xmlattr>.value")
                        +
                        ".parent",
                        path
                    );
                    __attribute_paths.put(
                        node->second.get<std::string>("<xmlattr>.value")
                        +
                        ".parent_type",
                        parent_type
                    );
                    __attribute_paths.put(
                        node->second.get<std::string>("<xmlattr>.value")
                        +
                        ".name",
                        element_name
                    );
                    if(parent_type==LINK || parent_type==ATTRIBUTE)
                        throw Error::Runtime(
                            "Attribute of a link or attribute "
                            "found in HDF5 file structure configuration!"
                        );
                    continue;
                } else if(type=="link") {
                    __link_paths.put(
                        node->second.get<std::string>("<xmlattr>.value"),
                        node_path
                    );
                    continue;
                } 
                fill_paths(node->second,
                           node_path,
                           get_node_type(node->second));
            }
        }
    }

    std::string SubPixHDF5File::substitute_ap_ind(
            const std::string &format_string,
            unsigned ap_ind)
    {
        size_t substitute_pos=format_string.find("%(ap_ind)");
        if(substitute_pos==std::string::npos) return format_string;
        char *formatted=new char[format_string.size()+100];
        std::string format_copy=format_string;
        format_copy.erase(substitute_pos+1, 8);
        int size=std::sprintf(formatted, format_copy.c_str(), ap_ind);
        if(size<0) throw Error::Runtime(
            "Failed to format aperture index name: '"
            +
            format_string
            +
            "'!"
        );
        std::string result(formatted, size);
        delete[] formatted;
        return result;
    }

    void SubPixHDF5File::configure(const char* xml_config_fname)
    {
#ifndef DEBUG
        H5::Exception::dontPrint();
#endif
        read_xml(xml_config_fname, __structure,
                 boost::property_tree::xml_parser::trim_whitespace);
        boost::property_tree::ptree::iterator root_node=__structure.begin();
#ifdef DEBUG
        boost::property_tree::ptree::iterator should_not_exist=root_node;
        assert(++should_not_exist==__structure.end());
#endif
        __structure=root_node->second;
        fill_paths();
    }

    H5::DataSet SubPixHDF5File::openDataSet(const std::string &dataset_key)
    {
        std::string path=__dataset_paths.get<std::string>(dataset_key, "");
        if(path=="") throw Error::InvalidArgument(
                "SubPixHDF5File::openDataSet",
                "No dataset is identified by '" + dataset_key + "'!"
        );
        try {
            return H5::H5File::openDataSet(path.c_str());
        } catch (H5::FileIException) {
            throw Error::HDF5NotFound(path, dataset_key, "dataset");
        } catch (H5::GroupIException) {
            throw Error::HDF5NotFound(path, dataset_key, "dataset");
        }
    }

    H5::Attribute SubPixHDF5File::openAttribute(
        const std::string &attribute_key
    )
    {
        boost::property_tree::ptree path=
            __attribute_paths.get_child(
                    attribute_key,
                    boost::property_tree::ptree()
            );
        if(path.size()==0) throw Error::InvalidArgument(
                "SubPixHDF5File::openAttribute",
                "No attribute is identified by '" + attribute_key + "'!"
        );
        __H5COMPONENT parent_type=static_cast<__H5COMPONENT>(
                path.get<int>("parent_type")
        );
        std::string parent_path=path.get<std::string>("parent");
        H5::H5Object *parent;
        H5::Group parent_group;
        H5::DataSet parent_dataset;
        switch(parent_type) {
            case GROUP:
                parent_group=H5::H5File::openGroup(parent_path);
                parent=&parent_group;
                break;
            case DATASET :
                parent_dataset=H5::H5File::openDataSet(parent_path);
                parent=&parent_dataset;
                break;
            default:  throw Error::Runtime("Invalid HDF5 file configuration "
                                           "encountered!");
        }
        std::string attribute_name=path.get<std::string>("name");
        if(!(parent->attrExists(attribute_name)))
            throw Error::HDF5NotFound(path.get<std::string>("parent")
                                      + "." + path.get<std::string>("name"),
                                      attribute_key,
                                      "attribute");
        return parent->openAttribute(path.get<std::string>("name"));
    }

    void SubPixHDF5File::create_attribute(H5::H5Object& destination,
                                          const std::string& name,
                                          const boost::any& value,
                                          int aperture_index)
    {
        if(value.empty()) return;
        const hsize_t single_value_dimensions[] = {1};
        const H5::DataSpace single_value_dataspace(1, 
                                                   single_value_dimensions);
        H5::DataType h5type;
        if(aperture_index<0 && value.type()==typeid(bool)) {
            H5::Attribute attribute=destination.createAttribute(
                    name,
                    H5::PredType::STD_U8LE,
                    single_value_dataspace
            );
            unsigned cast_value=(boost::any_cast<bool>(value) ? 1 : 0);
            attribute.write(H5::PredType::NATIVE_UINT, &cast_value);
        } else if(aperture_index<0 && value.type()==typeid(unsigned)) {
            H5::Attribute attribute=destination.createAttribute(
                    name,
                    H5::PredType::STD_U32LE,
                    single_value_dataspace
            );
            unsigned cast_value=boost::any_cast<unsigned>(value);
            attribute.write(H5::PredType::NATIVE_UINT, &cast_value);
        } else if(aperture_index<0 && value.type()==typeid(int)) {
            H5::Attribute attribute=destination.createAttribute(
                    name,
                    H5::PredType::STD_I32LE,
                    single_value_dataspace
            );
            int cast_value=boost::any_cast<int>(value);
            attribute.write(H5::PredType::NATIVE_INT, &cast_value);
        } else if(aperture_index<0 && value.type()==typeid(double)) {
            H5::Attribute attribute=destination.createAttribute(
                    name,
                    H5::PredType::IEEE_F32LE,
                    single_value_dataspace
            );
            double cast_value=boost::any_cast<double>(value);
            attribute.write(H5::PredType::NATIVE_DOUBLE, &cast_value);
        } else if(aperture_index<0 
                  &&
                  value.type()==typeid(std::pair<double, double>)) {
            std::pair<double, double> cast_value =
                boost::any_cast< std::pair<double, double> >(value);
            hsize_t size=2;
            std::valarray<double> data(2);
            data[0]=cast_value.first;
            data[1]=cast_value.second;
            H5::Attribute attribute=destination.createAttribute(
                    name,
                    H5::PredType::IEEE_F32LE,
                    H5::DataSpace(1, &size)
            );
            attribute.write(H5::PredType::NATIVE_DOUBLE,
                            &(data[0]));
        } else {
            try {
                OutputArray<unsigned> array(value);
                hsize_t size=(aperture_index<0 ? array.size() : 1);
                H5::Attribute attribute=destination.createAttribute(
                        name,
                        H5::PredType::STD_U32LE,
                        H5::DataSpace(1, &size)
                );
                attribute.write(H5::PredType::NATIVE_UINT,
                                array.data() + std::max(0, aperture_index));
                return;
            } catch(boost::bad_any_cast) {}
            try {
                OutputArray<int> array(value);
                hsize_t size=(aperture_index<0 ? array.size() : 1);
                H5::Attribute attribute=destination.createAttribute(
                        name,
                        H5::PredType::STD_I32LE,
                        H5::DataSpace(1, &size)
                );
                attribute.write(H5::PredType::NATIVE_INT,
                                array.data() + std::max(0, aperture_index));
                return;
            } catch(boost::bad_any_cast) {}
            try {
                OutputArray<double> array(value);
                hsize_t size=(aperture_index<0 ? array.size() : 1);
                H5::Attribute attribute=destination.createAttribute(
                        name,
                        H5::PredType::IEEE_F32LE,
                        H5::DataSpace(1, &size)
                );
                attribute.write(H5::PredType::NATIVE_DOUBLE,
                                array.data() + std::max(0, aperture_index));
                return;
            } catch(boost::bad_any_cast) {}
            if(aperture_index>0) throw Error::Type(
                    "Failed to convert the per-aperture attribute '"
                    + name
                    + "' to an array!"
            );
            try {
                std::string text=boost::any_cast<std::string>(value);
                H5::StrType AttributeType(H5::PredType::C_S1, text.size());
                H5::Attribute attribute= destination.createAttribute(
                        name,
                        AttributeType,
                        single_value_dataspace
                        );
                attribute.write(AttributeType, text.c_str());
            } catch(boost::bad_any_cast) {
                throw Error::Type("Failed to convert the attribute '"
                                  + name
                                  +	"' back to any useable type!");
            }
        }
    }

    void SubPixHDF5File::add_attributes(
            const boost::property_tree::ptree& structure_node,
            const H5IODataTree& data,
            H5::H5Object& destination,
            bool overwrite,
            int aperture_index)
    {
        using boost::property_tree::ptree;
        for(ptree::const_iterator sub_node=structure_node.begin();
                sub_node!=structure_node.end(); ++sub_node)
            if(		sub_node->first!="<xmlattr>"
                    &&
                    get_node_type(sub_node->second)==ATTRIBUTE) {
                std::string attribute_id=sub_node->second.get<std::string>(
                        "<xmlattr>.value"
                );

                if(attribute_id == "psffit.varnames") continue;
                const boost::any &attribute_value=
                    data.get<boost::any>(attribute_id, boost::any());
                if(attribute_value.empty()) continue;

                std::string attribute_name=sub_node->second.get(
                        "<xmlattr>.element_name",
                        sub_node->first
                );

                if(destination.attrExists(attribute_name)) {
                    if(overwrite) destination.removeAttr(attribute_name);
                    else throw Error::HDF5(
                            "",
                            "Attribute named '"
                            + attribute_name
                            + "' already exists!"
                    );
                }
                create_attribute(
                        destination,
                        attribute_name,
                        attribute_value,
                        aperture_index
                );
            }
    }

    H5::DSetCreatPropList SubPixHDF5File::compression_proplist(
        const std::string &compression_string,
        bool floating_point_data,
        hsize_t dataset_size)
    {
        if(compression_string.size()==0)
            return H5::DSetCreatPropList::DEFAULT;
        H5::DSetCreatPropList result;
        result.setChunk(1, (dataset_size<__chunk_size
                            ? &dataset_size
                            : &__chunk_size));
        std::istringstream compression_stream(compression_string);
        while(!compression_stream.eof()) {
            std::string to_apply;
            std::getline(compression_stream, to_apply, ';');
            if(to_apply=="shuffle") result.setShuffle();
            else {
                std::istringstream to_apply_stream(to_apply);
                std::string filter;
                std::getline(to_apply_stream, filter, ':');
                unsigned parameter;
                to_apply_stream >> parameter;
                if(filter=="gzip") result.setDeflate(parameter);
                else if(filter=="scaleoffset") {
                    std::valarray<unsigned> filter_arguments(2);
                    filter_arguments[0]=(floating_point_data
                                         ? H5Z_SO_FLOAT_DSCALE
                                         : H5Z_SO_INT);
                    filter_arguments[1]=parameter;
                    result.setFilter(H5Z_FILTER_SCALEOFFSET,
                                     H5Z_FLAG_MANDATORY,
                                     2,
                                     &(filter_arguments[0]));
                } else throw Error::InvalidArgument(
                        "H5IO.cpp: compression_proplist",
                        "Unrecognized compression filter: "+filter
                );
            }
        }
        return result;
    }

    H5::DataSet SubPixHDF5File::add_psfmap(const std::string& dataset_name,
                                           H5::Group& destination,
                                           const H5IODataTree& data,
                                           const std::string &compression,
                                           bool overwrite)
    {
        std::string psfmodel = data.get<std::string>("psffit.model",
                                                     "",
                                                     translate_string);
        assert(psfmodel != "");
        const boost::any&
            coefficients_any = data.get<boost::any>("psffit.psfmap",
                                                    boost::any());
        if(coefficients_any.empty()) return H5::DataSet();
        OutputArray<double> coefficients_array(coefficients_any);
        std::valarray<hsize_t> output_dimensions;
        if(psfmodel == "sdk") {
            assert(coefficients_array.size() % 3 == 0);
            output_dimensions.resize(2);
            output_dimensions[0] = 3;
            output_dimensions[1] = coefficients_array.size() / 3;
        } else if(psfmodel == "bicubic" || psfmodel == "zero") {
            std::string grid_str = data.get<std::string>("psffit.grid",
                                                         "",
                                                         translate_string);
            assert(grid_str != "");
            unsigned grid_split = grid_str.find_first_of(';');
            std::string grid_x_str = grid_str.substr(0, grid_split),
                        grid_y_str = grid_str.substr(grid_split+1);
            unsigned x_resolution = std::count(grid_x_str.begin(),
                                               grid_x_str.end(),
                                               ','),
                     y_resolution = std::count(grid_y_str.begin(),
                                               grid_y_str.end(),
                                               ',');
            assert(x_resolution >= 1);
            assert(y_resolution >= 1);
            output_dimensions.resize(4);
            output_dimensions[0] = 4;
            output_dimensions[1] = y_resolution - 1;
            output_dimensions[2] = x_resolution - 1;
            hsize_t num_psf_params = output_dimensions[0]
                                     *
                                     output_dimensions[1]
                                     *
                                     output_dimensions[2];
            assert(
                    (coefficients_array.size() == 0 && num_psf_params == 0)
                    ||
                    coefficients_array.size() % num_psf_params == 0
            );
            if (coefficients_array.size() == 0) output_dimensions[3] = 0;
            else output_dimensions[3] =
                coefficients_array.size() / num_psf_params;
        } else throw Error::InvalidArgument(
                "SubPixHDF5File::add_psfmap"
                ,
                "Unrecognized PSF model: " + psfmodel + " encountered while "
                "outputting the PSF map."
        );
        H5::DataSpace output_dataspace(output_dimensions.size(),
                                       &(output_dimensions[0]));
        H5::DSetCreatPropList creation_properties = compression_proplist(
                compression,
                true,
                coefficients_array.size()
        );

        if(overwrite) {
            try {destination.unlink(dataset_name);}
            catch(...) {}
        }
        H5::DataSet dataset = destination.createDataSet(
                dataset_name,
                dataset_type("psffit.psfmap"),
                output_dataspace,
                creation_properties
        );
        dataset.write(coefficients_array.data(), H5::PredType::NATIVE_DOUBLE);
        return dataset;
    }

    void SubPixHDF5File::add_psfmap_variable_names(
            const PSF::MapVarListType &variables
    )
    {
        boost::property_tree::ptree varname_path=__attribute_paths.get_child(
                "psffit.varnames",
                boost::property_tree::ptree()
        );
        if(varname_path.size()==0) return;
        __H5COMPONENT parent_type=static_cast<__H5COMPONENT>(
                varname_path.get<int>("parent_type")
        );
        std::string parent_path=varname_path.get<std::string>("parent");
        H5::H5Object *parent;
        H5::Group parent_group;
        H5::DataSet parent_dataset;
        switch(parent_type) {
            case GROUP:
                parent_group=H5::H5File::openGroup(parent_path);
                parent=&parent_group;
                break;
            case DATASET :
                parent_dataset=H5::H5File::openDataSet(parent_path);
                parent=&parent_dataset;
                break;
            default:  throw Error::Runtime("Invalid HDF5 file configuration "
                                           "encountered!");
        }
        std::string attribute_name=varname_path.get<std::string>("name");

        std::vector<const char*> variable_names(variables.size());
        std::vector<const char*>::iterator
            var_name_i = variable_names.begin();
        for(
                PSF::MapVarListType::const_iterator var_i = variables.begin();
                var_i != variables.end();
                ++var_i
        ) {
            *var_name_i = var_i->first.c_str();
            ++var_name_i;
        }
        hsize_t num_variables = variables.size();
        H5::Attribute attribute = parent->createAttribute(
                attribute_name,
                VARLEN_STR_TYPE,
                H5::DataSpace(1, &num_variables)
        );
        attribute.write(VARLEN_STR_TYPE, &(variable_names[0]));
    }

    H5::DataSet SubPixHDF5File::add_psfmap_variables(
                const std::string& dataset_name,
                H5::Group& destination,
                const H5IODataTree& data,
                const std::string &compression,
                bool overwrite
    )
    {
        IOTreeBase::const_assoc_iterator dataset_node = find_dataset_data(
                "psffit.variables", -1, data
        );
        if(dataset_node == data.not_found()) return H5::DataSet();
        const PSF::MapVarListType &variables = 
            TranslateToAny<PSF::MapVarListType>().get_value(
                dataset_node->second.data()
            );
        hsize_t num_vars = variables.size(),
                num_sources = variables.front().second.size();
        std::valarray<double> output_data(num_vars * num_sources);
        unsigned slice_start = 0;
        for(
                PSF::MapVarListType::const_iterator var_i = variables.begin();
                var_i != variables.end();
                ++var_i
        ) {
            assert(var_i->second.size() == num_sources);
            output_data[std::slice(slice_start, num_sources, 1)] = 
                var_i->second;
            slice_start += num_sources;
        }

        hsize_t output_dimensions[2] = {num_vars, num_sources};
        H5::DSetCreatPropList creation_properties = compression_proplist(
                compression,
                true,
                output_data.size()
        );
        if(overwrite) {
            try {destination.unlink(dataset_name);}
            catch(...) {}
        }
        H5::DataSet dataset = destination.createDataSet(
                dataset_name,
                dataset_type("psffit.variables"),
                H5::DataSpace(2, output_dimensions),
                creation_properties
        );
        dataset.write(&(output_data[0]), H5::PredType::NATIVE_DOUBLE);

        add_psfmap_variable_names(variables);

        return dataset;
    }

    IOTreeBase::const_assoc_iterator SubPixHDF5File::find_dataset_data(
            const std::string &dataset_id,
            int aperture_index,
            const H5IODataTree& data
    )
    {
        typedef IOTreeBase::path_type path;
        path dataset_path(dataset_id);
        if(aperture_index>=0) 
            dataset_path /= path(
                    static_cast<std::ostringstream*>( 
                        &(std::ostringstream() << aperture_index)
                    )->str()
            );
#ifdef DEBUG
        std::string dset_path_str = dataset_path.dump();
#endif
        IOTreeBase::const_assoc_iterator
            dataset_node = data.find(dataset_path.reduce()),
            not_found = data.not_found();
        while(!dataset_path.empty() && dataset_node != not_found) {
            not_found = dataset_node->second.not_found();
#ifdef DEBUG
            std::string key = dataset_path.reduce();
            std::cerr << "Getting sub-node '" << key << "' of '" 
                      << dataset_node->first << "'" << std::endl;
            dataset_node = dataset_node->second.find(key);
#else
            dataset_node = dataset_node->second.find(dataset_path.reduce());
#endif
        }
        if(
                dataset_node != not_found
                &&
                dataset_node->second.data().empty()
        ) {
#ifdef DEBUG
            std::cerr << "Checking for filename " << getFileName()
                      << " split under node " << dset_path_str << std::endl;
#endif
            not_found = dataset_node->second.not_found();
            dataset_node = dataset_node->second.find(getFileName());
        }
        if(
                dataset_node == not_found
                ||
                dataset_node->second.data().empty()
        ) {

#ifdef DEBUG
            std::cerr << "Value for dataset: " << dataset_id
                      << " not found in output data tree. Node "
                      << dset_path_str << " found: " 
                      << (dataset_node != data.not_found())
                      << std::endl;
#endif
            return data.not_found();
        } else {
#ifdef DEBUG
            std::cerr << "Found dataset: " << dataset_id << std::endl;
#endif
            return dataset_node;
        }
    }

    H5::DataSet SubPixHDF5File::add_generic_dataset(
            const std::string& dataset_name,
            H5::Group& destination,
            const H5::DataType &output_data_type,
            const boost::any &any_array,
            const boost::property_tree::ptree& structure_node,
            bool overwrite
    )
    {
        if(
                output_data_type == H5::PredType::STD_U8LE
                || output_data_type == H5::PredType::STD_U16LE
                || output_data_type == H5::PredType::STD_U32LE
        ) {
            return add_1d_dataset<unsigned>(dataset_name,
                                            destination,
                                            H5::PredType::NATIVE_UINT,
                                            output_data_type,
                                            false,
                                            any_array,
                                            structure_node,
                                            overwrite);
        } else if(
                output_data_type == H5::PredType::STD_I8LE
                || output_data_type == H5::PredType::STD_I16LE
                || output_data_type == H5::PredType::STD_I32LE
        ) {
            return add_1d_dataset<int>(dataset_name,
                                       destination,
                                       H5::PredType::NATIVE_INT,
                                       output_data_type,
                                       false,
                                       any_array,
                                       structure_node,
                                       overwrite);
        } else if(
                output_data_type == H5::PredType::IEEE_F32LE
                || output_data_type == H5::PredType::IEEE_F64LE
        ) {
            return add_1d_dataset<double>(dataset_name,
                                          destination,
                                          H5::PredType::NATIVE_DOUBLE,
                                          output_data_type,
                                          true,
                                          any_array,
                                          structure_node,
                                          overwrite);
        } else if(output_data_type == VARLEN_STR_TYPE) {
            return add_1d_dataset<char*>(dataset_name,
                                         destination,
                                         VARLEN_STR_TYPE,
                                         output_data_type,
                                         false,
                                         any_array,
                                         structure_node,
                                         overwrite);
        }
        throw Error::HDF5(
            "Unexpected data type found when building HDF5 file!"
        );
    }

    void SubPixHDF5File::add_dataset(
            const std::string& dataset_name,
            H5::Group& destination,
            const H5IODataTree& data,
            const boost::property_tree::ptree& structure_node,
            bool overwrite,
            int aperture_index
    )
    {
        std::string dataset_id = structure_node.get<std::string>(
                "<xmlattr>.value"
        );
        H5::DataSet dataset;
        if(
                dataset_id == "psffit.psfmap"
                &&
                data.get_optional<boost::any>("psffit.model")
        ) {
            assert(aperture_index == -1);
            dataset = add_psfmap(
                    dataset_name,
                    destination,
                    data,
                    structure_node.get<std::string>("<xmlattr>.compression",
                                                    ""),
                    overwrite
            );
        } else if(dataset_id == "psffit.variables") {
            assert(aperture_index == -1);
            dataset = add_psfmap_variables(
                    dataset_name,
                    destination,
                    data,
                    structure_node.get<std::string>("<xmlattr>.compression",
                                                    ""),
                    overwrite
            );
        } else {
            IOTreeBase::const_assoc_iterator dataset_node = 
                find_dataset_data(dataset_id, aperture_index, data);
            if(dataset_node != data.not_found()) {
                H5::DataType output_data_type = dataset_type(dataset_id);
                dataset = add_generic_dataset(dataset_name,
                                              destination,
                                              output_data_type,
                                              dataset_node->second.data(),
                                              structure_node,
                                              overwrite);
            }
        }
        add_attributes(structure_node,
                       data,
                       dataset,
                       overwrite,
                       aperture_index);
    }

    void SubPixHDF5File::add_link(const std::string &dataset_to_link_id,
                                  const std::string &link_name,
                                  H5::Group &destination,
                                  bool overwrite)
    {
        std::string target_path=__dataset_paths.get<std::string>(
                dataset_to_link_id
        );
        try {
            H5G_stat_t info;
            destination.getObjinfo(link_name.c_str(), false, info);
            if(info.type==H5G_LINK) {
                if(overwrite) destination.unlink(link_name);
                else throw Error::HDF5(
                    "A soft link named '" + link_name
                    + "' already exists, not overwriting!"
                );
            } else throw Error::HDF5(
                    "",
                    "On object named '" + link_name
                    + "' already exists and is not a soft link!");
        } catch(H5::GroupIException) {}
#ifndef NDEBUG
        std::cerr << "Creating link: " << link_name
                  << " -> "
                  << target_path
                  << " (" << dataset_to_link_id << ")"
                  << std::endl;
#endif
        destination.link(H5L_TYPE_SOFT, target_path, link_name);
    }

    void SubPixHDF5File::write_node(const H5IODataTree& data,
                                    H5::Group& destination,
                                    const std::string &new_name,
                                    const boost::property_tree::ptree& node,
                                    bool overwrite,
                                    int aperture_index)
    {
        __H5COMPONENT node_type=get_node_type(node);
        if(node_type==GROUP) {
            H5::Group group;
            try {
                group=destination.openGroup(new_name);
            } catch( H5::GroupIException ) {
                group=destination.createGroup(new_name);
            } catch ( H5::FileIException ) {
                group=destination.createGroup(new_name);
            }
            write_data(data, group, node, overwrite, aperture_index);
        } else if(node_type==DATASET) {
            add_dataset(new_name,
                        destination,
                        data,
                        node,
                        overwrite,
                        aperture_index);
            if(
                (
                    node.get<std::string>("<xmlattr>.value")
                    ==
                    "projsrc.srcid.name"
                )
                &&
                data.get<boost::any>(
                    "projsrc.srcid.name",
                    boost::any()
                ).empty()
                &&
                !data.get<boost::any>(
                    "projsrc.srcid.source", 
                    boost::any()
                ).empty()
            ) {
                H5::Group group;
                try {
                    group=destination.openGroup(new_name);
                } catch( H5::GroupIException ) {
                    group=destination.createGroup(new_name);
                } catch( H5::FileIException ) {
                    group=destination.createGroup(new_name);
                };
                write_data(data, group, node, overwrite, aperture_index);
            }
        } else if(node_type==LINK) {
            boost::any target=data.get<boost::any>(
                    node.get<std::string>("<xmlattr>.value"),
                    boost::any()
            );
            if(!target.empty()) add_link(
                    boost::any_cast<std::string>(target),
                    new_name,
                    destination,
                    overwrite
            );
        }
    }

    void SubPixHDF5File::write_data(
        const H5IODataTree &data,
        H5::Group &destination,
        const boost::property_tree::ptree &structure,
        bool overwrite,
        int aperture_index
    )
    {
        using boost::property_tree::ptree;
        add_attributes(structure,
                       data,
                       destination,
                       overwrite,
                       aperture_index);
        int max_aperture;
        if(aperture_index>=0) max_aperture=aperture_index;
        else {
            boost::any apertures=data.get<boost::any>("apphot.aperture",
                                                      boost::any());
            max_aperture=(apertures.empty()
                          ? -1 
                          : OutputArray<double>(apertures).size()-1);
        }
        for(ptree::const_iterator node=structure.begin();
                node!=structure.end(); ++node) {
            if(node->first!="<xmlattr>") {
                std::string element_name=node->second.get(
                        "<xmlattr>.element_name",
                        node->first
                );
                try {
                    int last_node_aperture=max_aperture;
                    for(
                            int current_ap_ind=aperture_index;
                            current_ap_ind<=last_node_aperture;
                            ++current_ap_ind
                    ) {
                        std::string new_name=substitute_ap_ind(
                                element_name,
                                std::max(0, current_ap_ind)
                        );
                        if(new_name==element_name)
                            last_node_aperture=aperture_index;
                        else current_ap_ind=std::max(current_ap_ind, 0);
                        write_node(data,
                                   destination,
                                   new_name,
                                   node->second,
                                   overwrite,
                                   current_ap_ind);
                    }
                } catch(Error::HDF5 &error) {
                    error.set_path(element_name+"/"+error.get_path());
                    throw;
                }
            }
        }
    }

    void SubPixHDF5File::read_pair_attribute(const H5::Attribute &attribute,
                                             H5IODataTree &data,
                                             const std::string &path)
    {
        H5::DataSpace data_space=attribute.getSpace();
#ifdef DEBUG
        hsize_t length;
        data_space.getSimpleExtentDims(&length);
        assert(length==2);
#endif
        H5T_class_t type_class=attribute.getTypeClass();
        if(type_class==H5T_INTEGER) {
            if(attribute.getIntType().getSign()==H5T_SGN_NONE) {
                unsigned values[2];
                attribute.read(H5::PredType::NATIVE_UINT, values);
                data.put(path,
                         std::pair<unsigned, unsigned>(values[0], values[1]),
                         TranslateToAny< std::pair<unsigned, unsigned> >());
            } else {
                int values[2];
                attribute.read(H5::PredType::NATIVE_INT, values);
                data.put(path,
                         std::pair<int, int>(values[0], values[1]),
                         TranslateToAny< std::pair<int, int> >());
            }
        } else if(type_class==H5T_FLOAT) {
            double values[2];
            attribute.read(H5::PredType::NATIVE_DOUBLE, values);
            data.put(path,
                     std::pair<double, double>(values[0], values[1]),
                     TranslateToAny< std::pair<double, double> >());
        } else throw Error::HDF5("Unrecognized data type for attribute "
                                 + path + " in " + getFileName());
    }

    void SubPixHDF5File::write(const H5IODataTree &data,
                               bool overwrite)
    {
        try {
            H5::Group root_group=openGroup("/");
            write_data(data, root_group, __structure, overwrite);
        } catch(Error::HDF5 &error) {
            error.set_path(getFileName()+"/"+error.get_path());
            throw;
        }
    }

    void SubPixHDF5File::read_psfmap(H5IODataTree &data)
    {
        std::vector<std::string> model(1, "psffit.model");
        read(model.begin(), model.end(), data);
        std::string path=__dataset_paths.get<std::string>("psffit.psfmap",
                                                          "");
        if(path=="") throw Error::Runtime("No configuration found for the "
                                          "path to the PSF map!");
        std::vector<std::string> terms(1, "psffit.terms");
        read(terms.begin(), terms.end(), data);

        if(
            data.get<std::string>("psffit.model", "", translate_string)
            ==
            "zero"
        ) {
            assert(!H5::H5File::nameExists(path));
            data.put("psffit.psfmap",
                     std::vector<double>(),
                     TranslateToAny< std::vector<double> >());
            return;
        }

        H5::DataSet coef_dataset = H5::H5File::openDataSet(path);
        assert(coef_dataset.getTypeClass() == H5T_FLOAT);
        H5::DataSpace file_space = coef_dataset.getSpace();
        int expected_dimensions;
#ifdef DEBUG
        unsigned expected_first_dimension;
#endif
        if(
                data.get<std::string>("psffit.model", "", translate_string)
                ==
                "sdk"
        ) {
            expected_dimensions = 2;
#ifdef DEBUG
            expected_first_dimension = 3;
#endif
        } else {
            expected_dimensions = 4;
#ifdef DEBUG
            expected_first_dimension = 4;
#endif
        }
        assert(file_space.getSimpleExtentNdims() == expected_dimensions);
        std::valarray<hsize_t> dimensions(expected_dimensions);
        file_space.getSimpleExtentDims(&(dimensions[0]));
#ifdef DEBUG
        assert(dimensions[0] == expected_first_dimension);
#endif
        file_space.selectAll();
        std::vector<double> result(file_space.getSimpleExtentNpoints());
        coef_dataset.read(&(result[0]), H5::PredType::NATIVE_DOUBLE);
        data.put("psffit.psfmap",
                 result,
                 TranslateToAny< std::vector<double> >());
    }

    void SubPixHDF5File::read_psfmap_variables(H5IODataTree &data)
    {
        H5::Attribute var_names_attribute = openAttribute("psffit.varnames");
        assert(var_names_attribute.getTypeClass() == H5T_STRING);
        hsize_t num_vars;
        var_names_attribute.getSpace().getSimpleExtentDims(&num_vars);
        std::vector<char*> varnames(num_vars);
        simple_read_from_h5(var_names_attribute,
                            var_names_attribute.getDataType(),
                            &(varnames[0]));

        H5::DataSet var_dataset = openDataSet("psffit.variables");
        assert(var_dataset.getTypeClass() == H5T_FLOAT);
        H5::DataSpace var_file_space = var_dataset.getSpace();
        assert(var_file_space.getSimpleExtentNdims() == 2);
        std::valarray<hsize_t> dimensions(2);
        var_file_space.getSimpleExtentDims(&(dimensions[0]));
        assert(dimensions[0] == num_vars);
        std::valarray<double> var_values(
                var_file_space.getSimpleExtentNpoints()
        );
        var_file_space.selectAll();
        var_dataset.read(&(var_values[0]), H5::PredType::NATIVE_DOUBLE);


        size_t num_sources = dimensions[1];
        PSF::MapVarListType variables;
        const double *variable_start = &(var_values[0]);
        for(unsigned var_index = 0; var_index < num_vars; ++var_index) {
            variables.push_back(
                    PSF::MapVariableType(
                        varnames[var_index],
                        std::valarray<double>(variable_start, num_sources)
                    )
            );
            variable_start += num_sources;
        }
        data.put("psffit.variables",
                 variables,
                 TranslateToAny<PSF::MapVarListType>());
    }

    /*int main(int argc, char** argv)
    {
        SubPixHDF5File::configure(argv[1]);
        H5IODataTree data_tree(argc, argv, "v0");
        TranslateToAny< std::valarray<double> > double_array_trans;
        TranslateToAny< std::valarray<unsigned> > unsigned_array_trans;
        std::valarray<double> x_array(121),
                              y_array(121),
                              apphot0(121),
                              apphot1(121),
                              apphot2(121),
                              apertures(3);
        apertures[0]=1.0;
        apertures[1]=2.5;
        apertures[2]=4.0;
        std::valarray<unsigned> field(121), source(121);
        unsigned index=0;
        for(double y=0; y<1.05; y+=0.1)
            for(double x=0; x<1.05; x+=0.1) {
                x_array[index]=x;
                y_array[index]=y;
                field[index]=static_cast<unsigned>(round(3*y));
                apphot0[index]=static_cast<unsigned>(y);
                apphot1[index]=static_cast<unsigned>(2*y);
                apphot2[index]=static_cast<unsigned>(4*y);
                source[index]=index;
                ++index;
            }
        data_tree.put("projsrc.x", &x_array, double_array_trans);
        data_tree.get< std::valarray<double> >("projsrc.x",
                                               std::valarray<double>(),
                                               double_array_trans);
        data_tree.put("projsrc.y", &y_array, double_array_trans);
        data_tree.put("projsrc.srcid.source", &source, unsigned_array_trans);
        data_tree.put("projsrc.srcid.field", &field, unsigned_array_trans);
        data_tree.put("projsrc.catalogue", "UCAC4", translate_string);
        data_tree.put("apphot.mag.0", apphot0, double_array_trans);
        data_tree.put("apphot.mag.1", apphot1, double_array_trans);
        data_tree.put("apphot.mag.2", apphot2, double_array_trans);
        data_tree.put("apphot.aperture", apertures, double_array_trans);
        std::ostringstream cmdline;
        for(int i=0; i<argc; ++i) {
            cmdline << "'" << argv[i] << "'";
            if(i!=argc-1) cmdline << " ";
        }
        data_tree.put("projsrc.cmdline", cmdline.str(), translate_string);
        SubPixHDF5File file("test.hdf5", H5F_ACC_TRUNC);
        file.write(data_tree, false);
        return 0;
    }*/

} //End IO namespace.
