#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <list>

#include "../FitPSF/CInterface.h"
#include "../IO/TranslateToAny.h"
#include "../IO/H5IODataTree.h"

///Create and set the configuration for PSF fitting.
FittingConfiguration *configure_psf_fitting(char *config_fname)
{
    std::ifstream config_file(config_fname);
    std::vector<std::string> config_keys(20), config_values(20);
    std::string equal;
    for(unsigned config_ind = 0; config_ind < 20; ++config_ind) {
        std::string config_line;
        std::getline(config_file, config_line);
        std::istringstream line_stream(config_line);
        line_stream >> config_keys[config_ind];
        line_stream >> equal;
        std::cerr << "Equal: " << equal << std::endl;
        assert(equal == "=");
        while(std::isspace(line_stream.get())) {}
        line_stream.unget();
        std::getline(line_stream, config_values[config_ind]);
    }

    FittingConfiguration *configuration = create_psffit_configuration();
    update_psffit_configuration(false, configuration,
                                "psf.model", "bicubic",
                                config_keys[0].c_str(), config_values[0].c_str(),
                                config_keys[1].c_str(), config_values[1].c_str(),
                                config_keys[2].c_str(), config_values[2].c_str(),
                                config_keys[3].c_str(), config_values[3].c_str(),
                                config_keys[4].c_str(), config_values[4].c_str(),
                                config_keys[5].c_str(), config_values[5].c_str(),
                                config_keys[6].c_str(), config_values[6].c_str(),
                                config_keys[7].c_str(), config_values[7].c_str(),
                                config_keys[8].c_str(), config_values[8].c_str(),
                                config_keys[9].c_str(), config_values[9].c_str(),
                                config_keys[10].c_str(), config_values[10].c_str(),
                                config_keys[11].c_str(), config_values[11].c_str(),
                                config_keys[12].c_str(), config_values[12].c_str(),
                                config_keys[13].c_str(), config_values[13].c_str(),
                                config_keys[14].c_str(), config_values[14].c_str(),
                                config_keys[15].c_str(), config_values[15].c_str(),
                                config_keys[16].c_str(), config_values[16].c_str(),
                                config_keys[17].c_str(), config_values[17].c_str(),
                                config_keys[18].c_str(), config_values[18].c_str(),
                                config_keys[19].c_str(), config_values[19].c_str(),
                                "");
    return configuration;
}

///\brief Read in the image pixel values from a file.
///
///The file should begin with two integer numbers specifying the x and y
///resolution of the image, followed by white space separated floating point
///values specifying the pixel values with the x coordinate varying faster.
double *read_image_data(char *filename,
                        unsigned long *x_resolution,
                        unsigned long *y_resolution)
{
    std::ifstream image_data_file(filename);
    image_data_file >> *x_resolution >> *y_resolution;
    double *image = new double[(*x_resolution) * (*y_resolution)];
    for(unsigned long y = 0; y < (*y_resolution); ++y)
        for(unsigned long x = 0; x < (*x_resolution); ++x) {
            double value;
            image_data_file >> value;
            image[x + (*x_resolution) * y] = value;
            assert(image_data_file);
        }
    return image;
}

///Print the pixels values of the given image to stdout.
template<typename DATA_TYPE>
void print_image(DATA_TYPE *image,
                 unsigned long x_resolution,
                 unsigned long y_resolution)
{
    std::cerr << "("
              << x_resolution << "x" << y_resolution
              << "):" << std::endl;

    for(unsigned long y = 0; y < y_resolution; ++y) {
        std::cerr << "    ";
        for(unsigned long x = 0; x < x_resolution; ++x)
            std::cerr << std::setw(10) << double(image[x + x_resolution * y]);
        std::cerr << std::endl;
    }

}

///\brief Create a mask with all good pixels.
char *create_clear_mask(unsigned long x_resolution, unsigned long y_resolution)
{
    char *mask = new char[x_resolution * y_resolution];
    for(unsigned long i = 0; i < x_resolution * y_resolution; ++i)
        mask[i] = 0.0;
    return mask;
}

///\brief Create an error image equal to sqrt(image).
double *create_error(double *image,
                     unsigned long x_resolution,
                     unsigned long y_resolution)
{
    double *error = new double[x_resolution * y_resolution];
    for(unsigned long i = 0; i < x_resolution * y_resolution; ++i)
        error[i] = std::sqrt(image[i]);
    return error;
}

///Allocate and fill a new C-style array with a list of strings.
void string_list_to_c_array(const std::list<std::string> &source,
                            char ***destination)
{
    *destination = new char*[source.size()];
    char **next_destination = *destination;
    for(
        std::list<std::string>::const_iterator src_i = source.begin();
        src_i != source.end();
        ++src_i
    ) {
        *next_destination = new char[src_i->size() + 1];
        std::strcpy(*next_destination, src_i->c_str());
        ++next_destination;
    }
}

///Allocate and fill a new C-style array with a list of doubles.
void double_list_to_c_array(const std::list<double> &source,
                            double *destination)
{
    for(
        std::list<double>::const_iterator src_i = source.begin();
        src_i != source.end();
        ++src_i
    ) {
        *destination = *src_i;
        ++destination;
    }
}

void read_column_names(std::istream &column_name_stream,
                       char ***column_names,
                       unsigned long *number_columns,
                       unsigned long *id_column_number)
{
    std::list<std::string> column_name_list;
    for(unsigned long column_number = 0; column_name_stream; ++column_number) {
        std::string name;
        column_name_stream >> name;
        if(column_name_stream) {
            if(name == "ID")
                *id_column_number = column_number;
            else
                column_name_list.push_back(name);
        }
    }
    *number_columns = column_name_list.size();
    string_list_to_c_array(column_name_list, column_names);
}

///\brief Read the list of sources to include in the fit.
void read_source_columns(char *filename,
                         char ***column_names,
                         char ***source_ids,
                         double **column_data,
                         unsigned long *number_columns,
                         unsigned long *number_sources)
{
    std::ifstream column_file(filename);
    std::string column_name_line;
    std::getline(column_file, column_name_line);
    std::istringstream column_name_stream(column_name_line);
    unsigned long id_column_number;
    read_column_names(column_name_stream,
                      column_names,
                      number_columns,
                      &id_column_number);
    std::vector< std::list<double> > column_vector(*number_columns);
    std::list<std::string> source_id_list;
    for(*number_sources = 0; column_file; ++(*number_sources)) {
        std::vector< std::list<double> >::iterator value_destination =
            column_vector.begin();
        for(
            unsigned long column_number = 0;
            column_number < (*number_columns) + 1;
            ++column_number
        ) {
            if(column_number == id_column_number) {
                std::string id;
                column_file >> id;
                if(column_file) {
                    std::cerr << "Read ID: " << id << std::endl;
                    source_id_list.push_back(id);
                }
            } else {
                double value;
                column_file >> value;
                if(column_file) {
                    std::cerr << "Read value " << value << std::endl;
                    value_destination->push_back(value);
                    value_destination++;
                }
            }
        }
    }
    string_list_to_c_array(source_id_list, source_ids);
    *number_sources = source_id_list.size();
    *column_data = new double[(*number_sources) * (*number_columns)];
    double *column_dest = *column_data;
    for(
        std::vector< std::list<double> >::iterator
            column_source = column_vector.begin();
        column_source != column_vector.end();
        ++column_source
    ) {
        double_list_to_c_array(*column_source,
                               column_dest);
        column_dest += *number_sources;
    }
}

BackgroundMeasureAnnulus *measure_backgrounds(unsigned long image_x_resolution,
                                              unsigned long image_y_resolution,
                                              double *image,
                                              double *error,
                                              char *mask,
                                              unsigned long number_columns,
                                              unsigned long number_sources,
                                              char **column_names,
                                              double *column_data)
{
    CoreImage *core_image = create_core_image(image_x_resolution,
                                              image_y_resolution,
                                              image,
                                              error,
                                              mask,
                                              true);
    BackgroundMeasureAnnulus *backgrounds = create_background_extractor(
        6.0,
        13.0,
        6.0,
        core_image,
        0.68
    );

    unsigned long x_column, y_column;

    for(unsigned long column_i = 0; column_i < number_columns; ++column_i) {
        if(strcmp(column_names[column_i], "x") == 0)
            x_column = column_i;
        else if(strcmp(column_names[column_i], "y") == 0)
            y_column = column_i;
    }
    add_source_list_to_background_extractor(
        backgrounds,
        column_data + x_column * number_sources,
        column_data + y_column * number_sources,
        number_sources
    );

    destroy_core_image(core_image);
    return backgrounds;
}

///Print the best fit coefficients in the given tree.
void print_coefficients(H5IODataTree *result_tree, unsigned long number_terms)
{
    double coefficients[4 * number_terms];
    query_result_tree(result_tree,
                      "psffit.psfmap.",
                      "[double]",
                      coefficients);
    for(unsigned param_i = 0; param_i < 4; ++param_i) {
        switch(param_i) {
            case 0: std::cerr << "f:"; break;
            case 1: std::cerr << "df_dx:"; break;
            case 2: std::cerr << "df_dy"; break;
            case 3: std::cerr << "d2f_dxdy"; break;
            default: assert(false);
        }
        std::cerr << std::endl;
        for(unsigned long term_i = 0; term_i < number_terms; ++term_i)
            std::cerr << coefficients[number_terms * param_i + term_i] << ", ";
        std::cerr << std::endl;
    }
}

int main(int argc, char **argv)
{
    unsigned long number_images = (argc - 2) / 2;
    unsigned long image_x_resolution,
                  image_y_resolution,
                  number_columns,
                  number_sources[number_images];
    double *images[number_images],
           *errors[number_images],
           *column_data[number_images],
           one = 1.0;
    char *masks[number_images];

    BackgroundMeasureAnnulus *backgrounds[number_images];

    char **column_names=NULL, **source_ids[number_images];

    std::cerr << "Processing " << number_images << " images." << std::endl;
    FittingConfiguration *configuration = configure_psf_fitting(
        argv[1]
    );

    for(unsigned image_index = 0; image_index < number_images; ++image_index) {
        images[image_index] = read_image_data(argv[2 + 2 * image_index],
                                              &image_x_resolution,
                                              &image_y_resolution);
        errors[image_index] = create_error(images[image_index],
                                           image_x_resolution,
                                           image_y_resolution);
        masks[image_index] = create_clear_mask(image_x_resolution,
                                               image_y_resolution);

        std::cerr << "Read image ";
        print_image(images[image_index],
                    image_x_resolution,
                    image_y_resolution);

        std::cerr << "Created mask ";
        print_image(masks[image_index],
                    image_x_resolution,
                    image_y_resolution);

        std::cerr << "Created error ";
        print_image(errors[image_index],
                    image_x_resolution,
                    image_y_resolution);

        char **image_column_names;
        unsigned long image_number_columns;
        read_source_columns(argv[3 + 2 * image_index],
                            &image_column_names,
                            &(source_ids[image_index]),
                            &(column_data[image_index]),
                            &image_number_columns,
                            &(number_sources[image_index]));
        if(column_names) {
            assert(image_number_columns == number_columns);
            for(unsigned col_i = 0; col_i < number_columns; ++col_i) {
                assert(
                    std::strcmp(column_names[col_i], image_column_names[col_i])
                    ==
                    0
                );
                delete[] image_column_names[col_i];
            }
            delete[] image_column_names;
        } else {
            number_columns = image_number_columns;
            column_names = image_column_names;
            std::cerr << "Column names (" << number_columns << "):" << std::endl;
            for(unsigned col_i = 0; col_i < number_columns; ++col_i)
                std::cerr << "    " << column_names[col_i] << std::endl;
        }

        std::cerr << "Sources ("
                  << number_sources[image_index]
                  << "):"
                  << std::endl;
        for(
            unsigned long src_i = 0;
            src_i < number_sources[image_index];
            ++src_i
        ) {
            std::cerr << "    " << source_ids[image_index][src_i];
            for(unsigned long col_i = 0; col_i < number_columns; ++col_i) {
                std::cerr << std::setw(25)
                          << column_data[
                                 image_index
                             ][
                                 src_i + col_i * number_sources[image_index]
                             ];
            }
            std::cerr << std::endl;
        }

        backgrounds[image_index] = measure_backgrounds(image_x_resolution,
                                                       image_y_resolution,
                                                       images[image_index],
                                                       errors[image_index],
                                                       masks[image_index],
                                                       number_columns,
                                                       number_sources[image_index],
                                                       column_names,
                                                       column_data[image_index]);
    }

    char zero_char=0;
    H5IODataTree *result_tree = create_result_tree(configuration, &zero_char);
    std::cerr << "Fit converged: "
              << piecewise_bicubic_fit(images,
                                       errors,
                                       masks,
                                       number_images,
                                       image_x_resolution,
                                       image_y_resolution,
                                       column_names,
                                       source_ids,
                                       column_data,
                                       number_sources,
                                       number_columns,
                                       backgrounds,
                                       configuration,
                                       &one,
                                       1,
                                       1,
                                       result_tree);

    print_coefficients(result_tree, 7);

    for(unsigned image_index = 0; image_index < number_images; ++image_index) {
        double *flux = new double[number_sources[image_index]];
        std::ostringstream image_index_str;
        image_index_str << image_index;
        query_result_tree(result_tree,
                          ("psffit.flux." + image_index_str.str()).c_str(),
                          "[double]",
                          flux);
        std::cerr << "Fluxes (sub-image: " << image_index << "):" << std::endl;
        for(
            unsigned long src_i = 0;
            src_i < number_sources[image_index];
            ++src_i
        )
            std::cerr << "    " << flux[src_i] << std::endl;

        const std::vector<std::string> &source_name_vector =
            reinterpret_cast<IO::H5IODataTree*>(
                result_tree
            )->get< std::vector<std::string> >(
                "projsrc.srcid.name." + image_index_str.str(),
                std::vector<std::string>(),
                IO::TranslateToAny< std::vector<std::string> >()
            );


        delete[] images[image_index];
        delete[] errors[image_index];
        delete[] masks[image_index];
        for(unsigned src_i = 0; src_i < number_sources[image_index]; ++src_i)
            delete[] source_ids[image_index][src_i];
        delete[] source_ids[image_index];
        delete[] column_data[image_index];
        destroy_background_extractor(backgrounds[image_index]);
        delete[] flux;
    }

    for(unsigned col_i = 0; col_i < number_columns; ++col_i)
        delete[] column_names[col_i];
    delete[] column_names;
    destroy_psffit_configuration(configuration);
    destroy_result_tree(result_tree);

    /*
    }
    return 0;*/
}
