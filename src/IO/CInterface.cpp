/**\file
 *
 * \brief The definitions of the functinos from CInterface.h
 *
 * \ingroup IO
 */

#define BUILDING_LIBRARY
#include "CInterface.h"
#include "H5IODataTree.h"
#include "OutputArray.h"
#include "../Core/Image.h"
#include <string>

const char MASK_OK = Core::MASK_OK;
const char MASK_CLEAR = Core::MASK_CLEAR;
const char MASK_FAULT = Core::MASK_FAULT;
const char MASK_HOT = Core::MASK_HOT;
const char MASK_COSMIC = Core::MASK_COSMIC;
const char MASK_OUTER = Core::MASK_OUTER;
const char MASK_OVERSATURATED = Core::MASK_OVERSATURATED;
const char MASK_LEAKED = Core::MASK_LEAKED;
const char MASK_SATURATED = Core::MASK_SATURATED;
const char MASK_INTERPOLATED = Core::MASK_INTERPOLATED;
const char MASK_BAD = Core::MASK_BAD;
const char MASK_ALL = Core::MASK_ALL;
const char MASK_NAN = Core::MASK_NAN;

void parse_hat_mask(const char *mask_string,
                    long x_resolution,
                    long y_resolution,
                    char *mask)
{
    return IO::parse_hat_mask(mask_string, x_resolution, y_resolution, mask);
}

H5IODataTree *create_result_tree(void *configuration, char *version_info)
{
    return reinterpret_cast<H5IODataTree*>(
        new IO::H5IODataTree(
            0,
            NULL,
            version_info,
            *reinterpret_cast<IO::CommandLineConfig*>(configuration)
        )
    );
}

void destroy_result_tree(H5IODataTree *tree)
{
#ifdef VERBOSE_DEBUG
    std::cerr << "Deleting H5IOData tree at " << tree << std::endl;
#endif
    delete reinterpret_cast<IO::H5IODataTree*>(tree);
#ifdef VERBOSE_DEBUG
    std::cerr << "Successfully Deleted H5IOData tree." << std::endl;
#endif
}

///Set result to a single value of type UNIT_TYPE.
template<typename UNIT_TYPE>
void get_single_value(
    ///The value to read from.
    const boost::any &value,

    ///The destination to write to.
    void *result
)
{
    *reinterpret_cast<UNIT_TYPE*>(result) = boost::any_cast<UNIT_TYPE>(value);
}

///Set result to a newly allocated c-style string.
void get_string_value(
    ///The value to read from.
    const boost::any &value,

    ///The destination to write to.
    void *result
)
{
    const std::string &source = IO::translate_string.get_value(value);
    char **destination = reinterpret_cast<char**>(result);
    *destination = reinterpret_cast<char*>(
        malloc(sizeof(char) * (source.size() + 1))
    );
    strcpy(*destination, source.c_str());
}

///Set result to a std::pair of values both of type UNIT_TYPE.
template<typename UNIT_TYPE>
void get_value_pair(
    ///The value to read from.
    const boost::any &value,

    ///The destination to write to.
    void *result
)
{
    UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);
    const std::pair<UNIT_TYPE, UNIT_TYPE> &source =
        IO::TranslateToAny< std::pair<UNIT_TYPE, UNIT_TYPE> >().get_value(value);

    destination[0] = source.first;
    destination[1] = source.second;
}

///Try copying a STL container of values, each of type UNIT_TYPE to result.
template<typename SOURCE_CONTAINER_TYPE, typename UNIT_TYPE>
bool try_copying_container(
    ///The value to read from.
    const boost::any &value,

    ///The destination to write to.
    void *result
)
{
    try {
        const SOURCE_CONTAINER_TYPE & input_container =
            IO::TranslateToAny<SOURCE_CONTAINER_TYPE>().get_value(value);

        UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);

        std::copy(input_container.begin(), input_container.end(), destination);

        return true;
    } catch(const boost::bad_any_cast &) {
        return false;
    }
}

///Try copying an array of values, each of type UNIT_TYPE to result.
template<typename SOURCE_ARRAY_TYPE, typename UNIT_TYPE>
bool try_copying_array(
    ///The value to read from.
    const boost::any &value,

    ///The destination to write to.
    void *result
)
{
    try {
        const SOURCE_ARRAY_TYPE &input_array =
            IO::TranslateToAny<SOURCE_ARRAY_TYPE>().get_value(value);

        UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);

        const UNIT_TYPE *start = &(input_array[0]),
                        *end = start + input_array.size();

        std::copy(start, end, destination);

        return true;
    } catch(const boost::bad_any_cast &) {
        return false;
    }
}

///Copy an array of values to a C-style array.
template<typename UNIT_TYPE> void copy_array(
    ///The input array of values to copy from.
    const boost::any &value,

    ///The destination to fill with the values. Must be pre-allocated.
    void *result
)
{
    UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);
    if(value.type() == typeid(UNIT_TYPE)) {
        *destination = boost::any_cast<const UNIT_TYPE&>(value);
        return;
    }
    typedef Eigen::Matrix<UNIT_TYPE, Eigen::Dynamic, 1> VectorEigen;
    typedef Eigen::Array<UNIT_TYPE, Eigen::Dynamic, 1> ArrayEigen;
#ifdef VERBOSE_DEBUG
    if(try_copying_container< std::vector<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                  result))
        return;
    else
        std::cerr << "Not vector" << std::endl;
    if(try_copying_container< std::list<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                result))
        return;
    else
        std::cerr << "Not list" << std::endl;
    if(try_copying_array< std::valarray<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                result))
        return;
    else
        std::cerr << "Not valarray" << std::endl;
    if(try_copying_array< VectorEigen, UNIT_TYPE >(value, result))
        return;
    else
        std::cerr << "Not Eigen::Vector" << std::endl;
    if(try_copying_array< ArrayEigen, UNIT_TYPE >(value, result))
        return;
    else
        std::cerr << "Not Eigen::Array" << std::endl;
    throw boost::bad_any_cast();
#else
    if(
        !(
            try_copying_container< std::vector<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                       result)
            ||
            try_copying_container< std::list<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                     result)
            ||
            try_copying_array< std::valarray<UNIT_TYPE>, UNIT_TYPE >(value,
                                                                     result)
            ||
            try_copying_array< VectorEigen, UNIT_TYPE >(value, result)
            ||
            try_copying_array< ArrayEigen, UNIT_TYPE >(value, result)
        )
    )
        throw boost::bad_any_cast();
#endif
}

///Copy an array of strings to a C-style array of char*.
void copy_string_array(
    ///The input array of values to copy from.
    const boost::any &value,

    ///The destination to fill with the values. Must be pre-allocated.
    void *result,

    ///The number of character in a single entry in result (entries are expected
    ///to be consecutive.
    int result_string_size
)
{
    IO::OutputArray<std::string> source_names(value);
    char *destination = reinterpret_cast<char*>(result);
    for(size_t i = 0; i < source_names.size(); ++i) {
#ifdef VERBOSE_DEBUG
        std::cerr << "Copying string " << source_names[i] << std::endl;
#endif
        strcpy(destination, source_names[i].c_str());
        destination += result_string_size;
    }
}

bool query_result_tree(H5IODataTree *tree,
                       const char *quantity,
                       const char *format,
                       void *result)
{

    const boost::any &value =
        reinterpret_cast<IO::H5IODataTree*>(tree)->get<boost::any>(quantity,
                                                                   boost::any());

    if(value.empty()) {
        std::cout << "Empty quantity: " << quantity << std::endl;
        return false;
    }

    if(strcmp(format, "str") == 0)
        get_string_value(value, result);
    else if(strcmp(format, "int") == 0)
        get_single_value<int>(value, result);
    else if(strcmp(format, "long") == 0)
        get_single_value<long>(value, result);
    else if(strcmp(format, "short") == 0)
        get_single_value<short>(value, result);
    else if(strcmp(format, "char") == 0)
        get_single_value<char>(value, result);
    else if(strcmp(format, "uint") == 0)
        get_single_value<unsigned>(value, result);
    else if(strcmp(format, "ulong") == 0)
        get_single_value<unsigned long>(value, result);
    else if(strcmp(format, "ushort") == 0)
        get_single_value<unsigned short>(value, result);
    else if(strcmp(format, "uchar") == 0)
        get_single_value<unsigned char>(value, result);
    else if(strcmp(format, "bool") == 0)
        get_single_value<bool>(value, result);
    else if(strcmp(format, "double") == 0)
        get_single_value<double>(value, result);
    else if(strcmp(format, "[int]") == 0)
        copy_array<int>(value, result);
    else if(strcmp(format, "[long]") == 0)
        copy_array<long>(value, result);
    else if(strcmp(format, "[short]") == 0)
        copy_array<short>(value, result);
    else if(strcmp(format, "[char]") == 0)
        copy_array<char>(value, result);
    else if(strcmp(format, "[uint]") == 0)
        copy_array<unsigned>(value, result);
    else if(strcmp(format, "[ulong]") == 0)
        copy_array<unsigned long>(value, result);
    else if(strcmp(format, "[ushort]") == 0)
        copy_array<unsigned short>(value, result);
    else if(strcmp(format, "[uchar]") == 0)
        copy_array<unsigned char>(value, result);
    else if(strcmp(format, "[bool]") == 0)
        copy_array<bool>(value, result);
    else if(strcmp(format, "[double]") == 0)
        copy_array<double>(value, result);
    else if(
        format[0] == '['
        &&
        format[1] == 'S'
    ) {
        char *last_character;
        long result_string_size = std::strtol(format + 2,
                                              &last_character,
                                              0);
        if(
            last_character != format + strlen(format) - 1
            ||
            *last_character != ']'
        ) {
            std::cerr << "Error: invalid argument to query_result_tree:"
                      << "invalid format: " << std::string(format)
                      << std::endl;
            return false;
        }
        copy_string_array(value, result, result_string_size);
    } else {
        int split_position = 0;
        while(format[split_position]!=':') {
            if(format[split_position] == '\0'){
                std::cerr << "Error: invalid argument to query_result_tree:"
                        << "invalid format: " << std::string(format)
                        << std::endl;
                return false;
            }
            ++split_position;
        }
        if(strncmp(format, format + split_position + 1, split_position) != 0){
            std::cerr << "Error: invalid argument to query_result_tree:"
                    << "invalid format: " << std::string(format)
                    << std::endl;
            return false;
        };
        if(strncmp(format, "int", split_position) == 0)
            get_value_pair<int>(value, result);
        else if(strncmp(format, "long", split_position) == 0)
            get_value_pair<long>(value, result);
        else if(strncmp(format, "short", split_position) == 0)
            get_value_pair<short>(value, result);
        else if(strncmp(format, "char", split_position) == 0)
            get_value_pair<char>(value, result);
        else if(strncmp(format, "uint", split_position) == 0)
            get_value_pair<unsigned>(value, result);
        else if(strncmp(format, "ulong", split_position) == 0)
            get_value_pair<unsigned long>(value, result);
        else if(strncmp(format, "ushort", split_position) == 0)
            get_value_pair<unsigned short>(value, result);
        else if(strncmp(format, "uchar", split_position) == 0)
            get_value_pair<unsigned char>(value, result);
        else if(strncmp(format, "bool", split_position) == 0)
            get_value_pair<bool>(value, result);
        else if(strncmp(format, "double", split_position) == 0)
            get_value_pair<double>(value, result);
        else {
            std::cerr << "Error: invalid argument to query_result_tree:"
                    << "invalid format: " << std::string(format)
                    << std::endl;
            return false;
        }
    }
    return true;
}

void update_result_tree(const char *quantity,
                        void *value,
                        const char *format,
                        unsigned length,
                        H5IODataTree *tree)
{
    reinterpret_cast<IO::H5IODataTree*>(tree)->add_c_array(quantity,
                                                           value,
                                                           format,
                                                           length);
}

LIB_PUBLIC bool get_psf_map_variables(H5IODataTree *output_data_tree,
                                      unsigned image_index,
                                      double *column_data)
{
    IO::H5IODataTree *real_output_data_tree =
        reinterpret_cast<IO::H5IODataTree*>(output_data_tree);

    std::ostringstream tree_path;
    tree_path << "psffit.variables." << image_index;
#ifdef VERBOSE_DEBUG
    std::cerr << "Getting PSF map variables from tree at "
              << real_output_data_tree
              << " with path "
              << tree_path.str()
              << std::endl;
#endif

    const PSF::MapVarListType &variables =
        real_output_data_tree->get<PSF::MapVarListType>(
            tree_path.str(),
            PSF::MapVarListType(),
            IO::TranslateToAny<PSF::MapVarListType>()
        );
    if(variables.size() == 0) {
#ifdef VERBOSE_DEBUG
        std::cerr << "Empty varibales entry in result tree!" << std::endl;
#endif
        return false;
    }
#ifdef VERBOSE_DEBUG
    else
        std::cerr << "Finished reading back psffit variables with size:"
                  << variables.size() << "x" << variables.begin()->second.size()
                  << std::endl;
#endif

    double *destination = column_data;
    for(
        PSF::MapVarListType::const_iterator var_i = variables.begin();
        var_i != variables.end();
        ++var_i
    ) {
        const double *start = &(var_i->second[0]),
                     *end = start + var_i->second.size();

        std::copy(start, end, destination);
        destination += var_i->second.size();
    }
    return true;
}

///\brief Add the names of all quantities that currently have a value in a
///data tree to a user supplied list.
void append_quantities_to_list(
    ///The tree to list the quantities of.
    const IO::IOTreeBase &tree,

    ///The list to add the names to.
    std::list<std::string> &quantities,

    ///The prefix to place before all entries (the path to the parent node).
    const std::string &prefix=""
)
{
    for(
        IO::IOTreeBase::const_iterator node_i = tree.begin();
        node_i != tree.end();
        ++node_i
    ) {
        if(!node_i->second.data().empty()) {
            quantities.push_back(prefix + node_i->first);
        }
        append_quantities_to_list(node_i->second,
                                  quantities,
                                  prefix + node_i->first + ".");
    }
}

void export_free(
    void *arg
)
{
    free(arg);
}

LIB_PUBLIC unsigned list_tree_quantities(H5IODataTree *tree,
                                         char ***quantities)
{
    IO::H5IODataTree *real_tree = reinterpret_cast<IO::H5IODataTree*>(tree);

    std::list<std::string> quantities_list;

    append_quantities_to_list(*real_tree, quantities_list);

    *quantities = new char*[quantities_list.size()];
    char **destination = *quantities;
    for(
        std::list<std::string>::const_iterator
            quantity_i = quantities_list.begin();
        quantity_i != quantities_list.end();
        ++quantity_i
    ) {
        *destination = reinterpret_cast<char*>(
            malloc(sizeof(char) * (quantity_i->size() + 1))
        );
        strcpy(*destination, quantity_i->c_str());
        ++destination;
    }
    return quantities_list.size();
}

LIB_PUBLIC void set_psf_map_variables(char **map_variable_names,
                                      double *map_variable_values,
                                      unsigned num_map_variables,
                                      unsigned num_sources,
                                      unsigned image_index,
                                      H5IODataTree *tree)
{
    PSF::MapVarListType variables;
    double *this_map_variable = map_variable_values;
    for(unsigned var_index = 0; var_index < num_map_variables; ++var_index) {
        variables.push_back(
            PSF::MapVariableType(
                map_variable_names[var_index],
                std::valarray<double>(this_map_variable, num_sources)
            )
        );
        this_map_variable += num_sources;
    }

    std::ostringstream tree_path;
    tree_path << "psffit.variables." << image_index;

    reinterpret_cast<IO::H5IODataTree*>(tree)->put<PSF::MapVarListType>(
        tree_path.str(),
        variables,
        IO::TranslateToAny<PSF::MapVarListType>()
    );
}
