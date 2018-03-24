#include "MapSourceContainer.h"

namespace PSF {

    MapSourceContainer::MapSourceContainer(const IO::H5IODataTree &data_tree,
                                           unsigned num_apertures)
    {
        const MapVarListType &variables =
            data_tree.get<MapVarListType>(
                "psffit.variables",
                MapVarListType(),
                IO::TranslateToAny<MapVarListType>()
            );

        std::vector<TermValarray> expansion_term_values;
        evaluate_term_expression(
                data_tree.get<std::string>("psffit.terms",
                                           "",
                                           IO::translate_string),
                variables.begin(),
                variables.end(),
                expansion_term_values
        );
        IO::OutputArray<double>
            x(data_tree.get<boost::any>("projsrc.x")),
            y(data_tree.get<boost::any>("projsrc.y")),
            amplitude(data_tree.get<boost::any>("psffit.amplitude")),
            background(data_tree.get<boost::any>("bg.value")),
            background_error(data_tree.get<boost::any>("bg.error"));
        IO::OutputArray<char*>
            source_name(data_tree.get<boost::any>("projsrc.srcid.name"));
        IO::OutputArray<unsigned>
            background_npix(data_tree.get<boost::any>("bg.npix"));

        unsigned num_sources = expansion_term_values[0].size(),
                 num_terms = expansion_term_values.size();
        if(data_tree.get<std::string>("psffit.model",
                                      "",
                                      IO::translate_string) == "zero")
            num_terms = 0;
        reserve(num_sources);

        for(unsigned src_index = 0; src_index < num_sources; ++src_index) {
            std::cerr << "source " << src_index
                      << ": " << source_name[src_index] << std::endl;
            push_back(
                MapSource(
                    Core::SourceID(source_name[src_index], true),
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
                new_terms[term_i] = expansion_term_values[term_i][src_index];
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
