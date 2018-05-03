/**\file
 *
 * \brief Utilities for manipulating HDF5 files.
 *
 * \ingroup IO
 */

#ifndef __HDF5_IO_H
#define __HDF5_IO_H

#include "../Core/SharedLibraryExportMacros.h"
#include "TranslateToAny.h"
#include "H5IODataTree.h"
#include "OutputArray.h"
#include "../PSF/Typedefs.h"
#include "../Core/Error.h"
#include "Eigen/Dense"
#include <H5Cpp.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/any.hpp>
#include <typeinfo>
#include <valarray>
#include <iostream>
#include <algorithm>

namespace IO {

    /* \brief A class for reading/writing HDF5 files imposing a user defined
     * structure.
     *
     * \ingroup IO
     */
    class LIB_PUBLIC SubPixHDF5File : public H5::H5File {
    private:
        ///Identifiers for the various HDF5 components.
        enum __H5COMPONENT{
            GROUP,		///< An HDF5 group
            DATASET,	///< An HDF5 dataset
            LINK,		///< An HDF5 link to a dataset
            ATTRIBUTE	///< An HDF5 attribute
        };

        ///\brief The type to use for storing variable length strings.
        ///
        ///Examples of quantities needing this type include:
        /// - named sources (not HAT-IDs).
        /// - PSF map variable names
        static const H5::StrType VARLEN_STR_TYPE;

        ///\brief The chunk size to use for chunked datasets
        ///
        ///All compressed datasets are chunked.
        static const hsize_t __chunk_size;

        ///The structure of the output file parsed on creation.
        static boost::property_tree::ptree __structure,
            
            ///The paths to the recognized datasets in __structure.
            __dataset_paths,
            
            ///The paths to the recognized attributes
            __attribute_paths,
            
            ///The paths to the various links
            __link_paths;

        ///The H5COMPONENT corresponding to a node in the HDF5 structure.
        static __H5COMPONENT get_node_type(
            ///The node whose type we need.
            const boost::property_tree::ptree &node
        );

        ///\brief Fills the __(dataset/property/link)_paths also fixed source 
        ///id entry.
        static void fill_paths(
            boost::property_tree::ptree& sub_structure=__structure,
            const std::string& parent_path="/",
            const __H5COMPONENT parent_type=GROUP
        );

        ///Returns the type to use for a dataset of the given name.
        const H5::DataType &dataset_type(
            ///The key identifying the dataset.
            const std::string &dataset_key
        );

        ///Substitutes the aperture index in a string.
        std::string substitute_ap_ind(
            ///The string to substitute in. It may or may not contain
            ///%(ap_ind)\<format\> style substitution.
            const std::string &format_string,

            ///The aperture index to assume.
            unsigned ap_ind
        );

        ///\brief Creates a new attribute, assuming no conflict with an 
        ///existing one.
        void create_attribute(
            ///The hdf5 location to add the attribute to.
            H5::H5Object& destination,

            ///The name of the new attribute to add.
            const std::string& name,

            ///The value of the attribute.
            const boost::any& value,

            ///The index of the aperture to output for aperture photometry
            ///dependent attributes.
            int aperture_index
        );

        ///Adds all known attributes to a group or a dataset.
        void add_attributes(
            ///The node in the structure corresponding to the location 
            ///being filled.
            const boost::property_tree::ptree& structure_node,

            ///The same as the data argument of SubPixHDF5File::write.
            const H5IODataTree& data,

            ///The location to which attributes must be added.
            H5::H5Object& destination,

            ///Should existing attributes be overwritten? If false,
            ///encountering an existintg attribute will throw an
            ///Error::HDF5 exception.
            bool overwrite,

            ///The index of the aperture to output for aperture photometry
            ///dependent attributes.
            int aperture_index
        );

        ///Adds the list of variables names attribute to the output file.
        void add_psfmap_variable_names(
            ///The list of variables being output.
            const std::list<PSF::MapVariableType> &variables
        );

        ///\brief Returns a H5::DSetCreatPropList that implements a given 
        ///compression.
        H5::DSetCreatPropList compression_proplist(
            ///What compression to use (if any) for a dataset. Can be a
            ///combination of: 'gzip:\<level\>', 'shuffle' and
            ///'scaleoffset:\<precision\>', separated by ';'. The precision for
            ///scale offset filters and the level of the gzip filter are
            ///defined the same as for the corresponding HDF5 filters.
            const std::string &compression_string,

            ///Will the dataset contain floating point (true) or integer
            ///(false) data
            bool floating_point_data,

            ///The size of the dataset being written (may determinene chunk
            ///size)
            hsize_t dataset_size
        );

        ///Creates (or overwwrites) the PSF map dataset.
        H5::DataSet add_psfmap(
            ///The name to give to the new dataset.
            const std::string& dataset_name,

            ///The group under which to place the dataset.
            H5::Group& destination,

            ///The same as the data argument to SubPixHDF5File::write
            const H5IODataTree& data,

            ///The compression string to use.
            const std::string &compression,

            ///If a dataset with the given name already exists, should it be
            ///overwritten? If false, an existing dataset throws Error::HDF5
            bool overwrite
        );

        ///\brief Creates (or overwwrites) the dataset containing the 
        ///variables required by the PSF map.
        H5::DataSet add_psfmap_variables(
            ///The name to give to the new dataset.
            const std::string& dataset_name,

            ///The group under which to place the dataset.
            H5::Group& destination,

            ///The same as the data argument to SubPixHDF5File::write
            const H5IODataTree& data,

            ///The compression string to use.
            const std::string &compression,

            ///If a dataset with the given name already exists, should it be
            ///overwritten? If false, an existing dataset throws Error::HDF5
            bool overwrite
        );

        ///\brief Return an iterator to the node containing the data for a
        ///dataset.
        ///
        ///Takes care of aperture index if necessary and split by output
        ///filename.
        IOTreeBase::const_assoc_iterator find_dataset_data(
            ///The ID of the dataset (e.g. projsrc.x).
            const std::string &dataset_id,

            ///See same name argument to add_dataset.
            int aperture_index,

            ///The output data tree to find the data in.
            const H5IODataTree& data
        );

        template<typename DATA_TYPE>
            H5::DataSet add_1d_dataset(
                const std::string& dataset_name,
                H5::Group& destination,
                const H5::DataType &memory_data_type,
                const H5::DataType &output_data_type,
                bool floating_point_data,
                const boost::any &any_array,
                const boost::property_tree::ptree& structure_node,
                bool overwrite
            );

        ///\brief Creates (or overwrites) a dataset other than the PSF map and
        ///PSF map variables.
        H5::DataSet add_generic_dataset(
            ///The name to give to the new dataset.
            const std::string& dataset_name,

            ///The group under which to place the dataset.
            H5::Group& destination,

            ///The data type to use when saving the dataset.
            const H5::DataType &output_data_type,

            ///The array of values to write.
            const boost::any &any_array,

            ///The node in the HDF5 file structure corresponding to the
            ///dataset.
            const boost::property_tree::ptree& structure_node,

            ///If a dataset with the given name already exists, should it be
            ///overwritten? If false, an existing dataset throws Error::HDF5
            bool overwrite
        );

        ///Creates (or overwrites) and fills a dataset.
        void add_dataset(
            ///The name to give to the new dataset.
            const std::string& dataset_name,

            ///The group under which to place the dataset.
            H5::Group& destination,

            ///The same as the data argument to SubPixHDF5File::write
            const H5IODataTree& data,

            ///The node in the HDF5 file structure corresponding to the
            ///dataset.
            const boost::property_tree::ptree& structure_node,

            ///If a dataset with the given name already exists, should it be
            ///overwritten? If false, an existing dataset throws Error::HDF5
            bool overwrite,

            ///The index of the aperture to output for aperture photometry
            ///dependent quantities.
            int aperture_index
        );

        ///Creates (or overwrites) a soft link to a dataset.
        void add_link(
            ///The ID of the dataset to link.
            const std::string &dataset_to_link_id,

            ///The name to give to the new link,
            const std::string &link_name,

            ///The group to create the link under.
            H5::Group &destination,

            ///Should a pre-existing link be overwritten? If false and a link
            ///with the given name already exists Error::HDF5 is thrown. If
            ///anything other than soft a link with the given name exists,
            ///Error::HDF5 is thrown.
            bool overwrite
        );

        ///Handles the writing of a single node from write_data.
        void write_node(
            ///See write_data.
            const H5IODataTree& data,

            ///See write_data.
            H5::Group& destination,

            ///The name to give to the newly created component.
            const std::string &new_name,

            ///How to structure the data one of the nodes under structure in
            ///write_data.
            const boost::property_tree::ptree& node,

            ///The same as the overwrite argument of SubPixHDF5File::write.
            bool overwrite,

            ///The index of the aperture to output for aperture photometry
            ///entries. Should be -1 if not currently outputting per aperture
            ///quantities.
            int aperture_index
        );

        ///Adds/overwrites data to (part of) the HDF5 file
        void write_data(
            ///The same as the data argument of SubPixHDF5File::write.
            const H5IODataTree& data,

            ///The location where the data should be written, should be a
            ///group in this file, possibly '/'.
            H5::Group& destination,

            ///How to structure the data (subtree of structure argument of
            ///SubPixHDF5File::write)>
            const boost::property_tree::ptree& structure,

            ///The same as the overwrite argument of SubPixHDF5File::write.
            bool overwrite,

            ///The index of the aperture to output for aperture photometry
            ///entries. Should be -1 if not currently outputting per aperture
            ///quantities.
            int aperture_index=-1
        );

        template<class DSET_ATTR, typename UNIT_TYPE>
        void read_1d(
            ///The dataset or attribute to read.
            const DSET_ATTR &source,

            ///The space of the dataset or attribute.
            const H5::DataSpace &data_space,

            ///The memory type to use when reading.
            const H5::DataType &memory_type,

            ///The data tree to fill.
            H5IODataTree &data,

            ///The path within the data tree to fill.
            const std::string &path
        );

        template<class DSET_ATTR>
        void read_1d_string(
            ///The dataset or attribute to read.
            const DSET_ATTR &source,

            ///The space of the dataset or attribute.
            const H5::DataSpace &data_space,

            ///The memory type to use when reading.
            const H5::DataType &memory_type,

            ///The data tree to fill.
            H5IODataTree &data,

            ///The path within the data tree to fill.
            const std::string &path
        );

        ///\brief Reading a scalar 1d dataset or attribute (return false if 
        ///not scalar).
        template<class DSET_ATTR>
        bool read_scalar(
            ///The dataset or attribute to read.
            const DSET_ATTR &source,

            ///The data tree to fill.
            H5IODataTree &data,

            ///The path within the data tree to fill.
            const std::string &path
        );

        ///Reads an attribute consisting of a pair of values.
        void read_pair_attribute(
            ///The attribute to read.
            const H5::Attribute &attribute,

            ///The data tree to fill
            H5IODataTree &data,

            ///The path within the data tree to fill.
            const std::string &path
        );
    public:
        ///Default constructor allows opening the file later.
        SubPixHDF5File() : H5File() {}

        ///Expose H5File constructor (see HDF5 documentation for description).
        SubPixHDF5File(const char *name,
                       unsigned int flags,
                       const H5::FileCreatPropList& create_plist =
                                H5::FileCreatPropList::DEFAULT,
                       const H5::FileAccPropList& access_plist =
                                H5::FileAccPropList::DEFAULT)
            : H5File(name, flags, create_plist, access_plist) {}

        ///Reads the structure of the desired HDF5 output from an XML file.
        static void configure(
            ///The name of the file to read the structure from.
            const char* xml_config_fname
        );

        ///Opens the dataset identified by the given key.
        H5::DataSet openDataSet(
            ///The key identifying the dataset (e.g. 'psffit.mag')
            const std::string &dataset_key
        );

        ///Opens the attribute identified by the given key.
        H5::Attribute openAttribute(
            ///The key identifying the attribute (e.g. 'psffit.cmdline')
            const std::string &attribute_key
        );

        ///Adds/overwrites data to an open file (attributes groups and links).
        void write(
            ///The content to write. Should be organized in a boost property
            ///tree with string keys and boost::any values. All keys should
            ///have one of the recognized names.
            const H5IODataTree &data,

            ///If an element in data is found to already exist this argument
            ///determines if it will be overwritten (True) or an
            ///Error::HDF5 exception will be raised (False).
            bool overwrite
        );

        ///Reads in the coefficients of the PSF map.
        void read_psfmap(
            ///The data tree to write the PSF map coefficients to.
            H5IODataTree &data
        );

        ///Reads in the variables which participate in the PSF map.
        void read_psfmap_variables(
            ///The data tree to write the PSF map variables to.
            H5IODataTree &data
        );

        ///Reads from the file into an I/O tree.
        template<class InputIterator>
        void read(
            ///The first of a set of string objects giving the keys to the
            ///quantities which should be read.
            InputIterator first,

            ///One past the last quantity to read.
            InputIterator last,

            ///The I/O tree to fill.
            H5IODataTree &data,

            ///Must all quantities be present. If false, missing quantities
            ///are simply not inserted in data.
            bool require_all = true
        );

    }; //End SubPixHDF5File class.

    ///Present a uniform reading function for attributes and datasets.
    LIB_PUBLIC void simple_read_from_h5(const H5::DataSet &dataset,
                                        const H5::DataType &memory_type,
                                        void *destination);

    ///Present a uniform reading function for attributes and datasets.
    LIB_PUBLIC void simple_read_from_h5(const H5::Attribute &attribute,
                                        const H5::DataType &memory_type,
                                        void *destination);

    template<typename DATA_TYPE>
        H5::DataSet SubPixHDF5File::add_1d_dataset(
            const std::string& dataset_name,
            H5::Group& destination,
            const H5::DataType &memory_data_type,
            const H5::DataType &output_data_type,
            bool floating_point_data,
            const boost::any &any_array,
            const boost::property_tree::ptree& structure_node,
            bool overwrite
        )
        {
            std::string compression = structure_node.get<std::string>(
                    "<xmlattr>.compression",
                    ""
            );
            OutputArray<DATA_TYPE> array(any_array);
#ifdef DEBUG
            std::cerr << "Writing array of size = "
                      << array.size() << std::endl
                      << "first 3 entries: "
                      << std::endl
                      << array[0] << std::endl
                      << array[1] << std::endl
                      << array[2] << std::endl;
#endif
            const void *array_data = array.data();
            H5::DataSpace dataspace = H5::DataSpace(1, &array.size());
            H5::DSetCreatPropList creation_properties = compression_proplist(
                compression,
                floating_point_data,
                array.size()
            );
#ifdef DEBUG
            std::cerr << "Creating " << dataset_name << " dataset with " 
                      << dataspace.getSimpleExtentNpoints() << " points!" 
                      << std::endl;
#endif
            if(overwrite) {
                try {destination.unlink(dataset_name);}
                catch(...) {}
            }
            H5::DataSet dataset = destination.createDataSet(
                dataset_name,
                output_data_type,
                dataspace,
                creation_properties
            );
            dataset.write(array_data, memory_data_type);
            return dataset;
        }

    template<class DSET_ATTR, typename UNIT_TYPE>
        void SubPixHDF5File::read_1d(const DSET_ATTR &source,
                                     const H5::DataSpace &data_space,
                                     const H5::DataType &memory_type,
                                     H5IODataTree &data,
                                     const std::string &path)
        {
            hsize_t length;
            data_space.getSimpleExtentDims(&length);
            std::vector<UNIT_TYPE> *result=new std::vector<UNIT_TYPE>(length);
            simple_read_from_h5(source, memory_type, &((*result)[0]));
            if(length>1)
                data.put(path,
                         result,
                         TranslateToAny< std::vector<UNIT_TYPE> >());
            else if(length==1)
                data.put(path,
                         (*result)[0],
                         TranslateToAny< UNIT_TYPE >());
            else
                throw Error::HDF5(getFileName(),
                                  "Zero length dataset encountered!");
        }

    template<class DSET_ATTR>
        void SubPixHDF5File::read_1d_string(const DSET_ATTR &source,
                                            const H5::DataSpace &data_space,
                                            const H5::DataType &memory_type,
                                            H5IODataTree &data,
                                            const std::string &path)
        {
            hsize_t length;
            data_space.getSimpleExtentDims(&length);
            if(length > 1) {
                std::vector<char*> result(length);
                simple_read_from_h5(source, memory_type, &(result[0]));
                data.put(path,
                         result,
                         TranslateToAny< std::vector<char*> >());
            } else if(length == 1) {
                H5::StrType datatype = source.getStrType();
                char *result = new char[datatype.getSize() + 1];
                simple_read_from_h5(source, memory_type, result);
                result[datatype.getSize()] = '\0';
                data.put(path,
                         std::string(result),
                         TranslateToAny<std::string>());
                delete[] result;
            } else 
                throw Error::HDF5(getFileName(),
                                  "Zero length dataset encountered!");
        }

    template<class DSET_ATTR>
        bool SubPixHDF5File::read_scalar(const DSET_ATTR &source,
                                         H5IODataTree &data,
                                         const std::string &path)
        {
            assert(source.getSpace().getSimpleExtentNdims() == 1);
            try {
                H5T_class_t type_class = source.getTypeClass();
                if(type_class == H5T_INTEGER) {
                    if(source.getIntType().getSign() == H5T_SGN_NONE)
                        read_1d<DSET_ATTR, unsigned>(
                            source,
                            source.getSpace(),
                            H5::PredType::NATIVE_UINT,
                            data,
                            path
                        );
                    else
                        read_1d<DSET_ATTR, int>(
                            source,
                            source.getSpace(),
                            H5::PredType::NATIVE_INT,
                            data,
                            path
                        );
                } else if(type_class==H5T_FLOAT) {
                    read_1d<DSET_ATTR, double>(
                        source,
                        source.getSpace(),
                        H5::PredType::NATIVE_DOUBLE,
                        data,
                        path
                    );								
                } else if(type_class==H5T_STRING) {
                    read_1d_string<DSET_ATTR>(
                        source,
                        source.getSpace(),
                        source.getDataType(),
                        data,
                        path
                    );
                } else
                    return false;
                return true;
            } catch(Error::HDF5 error) {
                error.set_path(error.get_path() + "/" + path);
                throw;
            }
        }

    template<class InputIterator>
        void SubPixHDF5File::read(InputIterator first,
                                  InputIterator last,
                                  H5IODataTree &data,
                                  bool require_all)

        {
            for(;first!=last; ++first) {
                if((*first) == "psffit.psfmap") {
                    read_psfmap(data);
                } else if ((*first) == "psffit.variables") {
                    read_psfmap_variables(data);
                } else {
                    std::string path = __dataset_paths.get<std::string>(
                        *first,
                        ""
                    );
                    try {
                        if(path=="") {
                            if(
                                !read_scalar(openAttribute(*first),
                                             data,
                                             *first)
                            )
                                throw Error::HDF5(
                                    "Attribute identified by '" + *first
                                    + "' in " + getFileName()
                                    + "' has an unrecognized data type!"
                                );
                        } else {
                            if(
                                !read_scalar(openDataSet(*first),
                                             data,
                                             *first)
                            )
                                throw Error::HDF5(
                                    "Dataset identified by '" + *first
                                    + "':  '" + getFileName() + "/" + path
                                    + "' has an unrecognized data type!"
                                );
                        }
                    } catch(Error::HDF5NotFound) {
                        if(require_all) throw;
                    }
                }
            }
        }

} //End IO namespace.

#endif
