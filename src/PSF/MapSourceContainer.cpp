#include "MapSourceContainer.h"

namespace PSF {

    MapSourceContainer::MapSourceContainer(
        const IO::H5IODataTree &data_tree,
        unsigned num_apertures,
        const std::string &data_tree_image_id
    )
    {
        std::string tree_suffix = (data_tree_image_id.empty()
                                   ? ""
                                   : "." + data_tree_image_id);

        std::string amplitude_name = (
            data_tree.count("psffit.amplitude" + tree_suffix)
            ? "psffit.amplitude"
            : "psffit.flux"
        ) + tree_suffix;
        IO::OutputArray<double>
            x(
                data_tree.get<boost::any>("projsrc.x" + tree_suffix)
            ),
            y(
                data_tree.get<boost::any>("projsrc.y" + tree_suffix)
            ),
            amplitude(
                data_tree.get<boost::any>(amplitude_name)
            ),
            background(
                data_tree.get<boost::any>("bg.value" + tree_suffix)
            ),
            background_error(
                data_tree.get<boost::any>("bg.error" + tree_suffix)
            ),
            expansion_term_values(
                data_tree.get<boost::any>("psffit.terms" + tree_suffix)
            );

        unsigned num_sources = x.size(),
                 num_terms = expansion_term_values.size() / num_sources;

        assert(y.size() == num_sources);
        assert(amplitude.size() == num_sources);
        assert(background.size() == num_sources);
        assert(background_error.size() == num_sources);
        assert(expansion_term_values.size() == num_sources * num_terms);

        IO::OutputArray<std::string> *source_name = NULL;
        IO::OutputArray<unsigned> *source_field = NULL,
                                  *source_id = NULL;
        if(data_tree.get_child("projsrc.srcid.name").count("0")) {
            source_name = new IO::OutputArray<std::string>(
                data_tree.get<boost::any>("projsrc.srcid.name" + tree_suffix)
            );
            assert(source_name->size() == num_sources);
        } else {
            source_field = new IO::OutputArray<unsigned>(
                data_tree.get<boost::any>("projsrc.srcid.field" + tree_suffix)
            );
            source_id = new IO::OutputArray<unsigned>(
                data_tree.get<boost::any>("projsrc.srcid.source" + tree_suffix)
            );

            assert(source_field->size() == num_sources);
            assert(source_id->size() == num_sources);
        }
        IO::OutputArray<unsigned>
            background_npix(
                data_tree.get<boost::any>("bg.npix" + tree_suffix)
            );
        assert(background_npix.size() == num_sources);


        assert(
            (
                data_tree.get<std::string>("psffit.model",
                                           "",
                                           IO::translate_string) != "zero"
            )
            ||
            num_terms == 0
        );
        reserve(num_sources);

        for(unsigned src_index = 0; src_index < num_sources; ++src_index) {
            push_back(
                MapSource(
                    (
                        source_name
                        ? Core::SourceID((*source_name)[src_index], true)
                        : Core::SourceID((*source_field)[src_index],
                                         (*source_id)[src_index])
                    ),
                    num_apertures,
                    x[src_index],
                    y[src_index],
                    Background::Source(
                        background[src_index] / amplitude[src_index],
                        background_error[src_index] / amplitude[src_index],
                        background_npix[src_index]
                    )
                )
            );
            Eigen::VectorXd &new_terms = back().expansion_terms();
            new_terms.resize(num_terms);
            for(unsigned term_i = 0; term_i < num_terms; ++term_i)
                new_terms[term_i] = expansion_term_values[src_index * num_terms
                                                          +
                                                          term_i];
        }
    }

    const std::set<std::string>&
        MapSourceContainer::required_data_tree_quantities()
    {
        const std::string additional_required_data_tree_quantities[] = {
            "projsrc.srcid.name",
            "projsrc.x",
            "projsrc.y",
            "bg.value",
            "bg.error",
            "bg.npix",
            "psffit.amplitude",
            "psffit.variables",
            "psffit.terms",
            "psffit.model"
        };

        static const std::set<std::string> required_quantities(
                    additional_required_data_tree_quantities,
                    additional_required_data_tree_quantities
                    +
                    sizeof(additional_required_data_tree_quantities)
                    /
                    sizeof(additional_required_data_tree_quantities[0])
        );
        return required_quantities;
    }

}
