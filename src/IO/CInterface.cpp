/**\file
 *
 * \brief The definitions of the functinos from CInterface.h
 *
 * \ingroup IO
 */

#define BUILDING_LIBRARY
#include "CInterface.h"
#include "H5IODataTree.h"
#include "../Core/Image.h"

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
    delete reinterpret_cast<IO::H5IODataTree*>(tree);
}

///Set result to a single value of type UNIT_TYPE.
template<typename UNIT_TYPE>
void get_single_value(const boost::any &value, void *result)
{
    *reinterpret_cast<UNIT_TYPE*>(result) = boost::any_cast<UNIT_TYPE>(value);
}

///Set result to a newly allocated c-style string.
void get_string_value(const boost::any &value, void *result)
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
void get_value_pair(const boost::any &value, void *result)
{
    UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);
    const std::pair<UNIT_TYPE, UNIT_TYPE> &source =
        IO::TranslateToAny< std::pair<UNIT_TYPE, UNIT_TYPE> >().get_value(value);

    destination[0] = source.first;
    destination[1] = source.second;
}

///Try copying a STL container of values, each of type UNIT_TYPE to result.
template<typename SOURCE_CONTAINER_TYPE, typename UNIT_TYPE>
bool try_copying_container(const boost::any &value, void *result)
{
    try {
        const SOURCE_CONTAINER_TYPE & input_container =
            IO::TranslateToAny<SOURCE_CONTAINER_TYPE>().get_value(value);

        UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);

        std::copy(input_container.begin(), input_container.end(), destination);

        return true;
    } catch(boost::bad_any_cast) {
        return false;
    }
}

///Try copying an array of values, each of type UNIT_TYPE to result.
template<typename SOURCE_ARRAY_TYPE, typename UNIT_TYPE>
bool try_copying_array(const boost::any &value, void *result)
{
    try {
        const SOURCE_ARRAY_TYPE &input_array =
            IO::TranslateToAny<SOURCE_ARRAY_TYPE>().get_value(value);

        UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);

        const UNIT_TYPE *start = &(input_array[0]),
                        *end = start + input_array.size();

        std::copy(start, end, destination);

        return true;
    } catch(boost::bad_any_cast) {
        return false;
    }
}

template<typename UNIT_TYPE> void copy_array(const boost::any &value, void *result)
{
    UNIT_TYPE *destination = reinterpret_cast<UNIT_TYPE*>(result);
    if(value.type() == typeid(UNIT_TYPE)) {
        *destination = boost::any_cast<const UNIT_TYPE&>(value);
        return;
    }
    typedef Eigen::Matrix<UNIT_TYPE, Eigen::Dynamic, 1> VectorEigen;
    typedef Eigen::Array<UNIT_TYPE, Eigen::Dynamic, 1> ArrayEigen;
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
    else {
        int split_position = 0;
        while(format[split_position]!=':') {
            if(format[split_position] == '\0')
                throw Error::InvalidArgument(
                    "query_result_tree",
                    "invalid format: " + std::string(format)
                );
            ++split_position;
        }
        if(strncmp(format, format + split_position + 1, split_position) != 0)
            throw Error::InvalidArgument(
                "query_result_tree",
                "invalid format: " + std::string(format)
            );
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
        else
            throw Error::InvalidArgument(
                "query_result_tree",
                "invalid format: " + std::string(format)
            );
    }
    return true;
}
