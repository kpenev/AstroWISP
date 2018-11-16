#ifndef __TERM_CALCULATOR_H
#define __TERM_CALCULATOR_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Typedefs.h"
#include "../Core/Error.h"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <valarray>
#include <sstream>

namespace PSF {

    ///\brief Same as valarray<double> but on assignment resizes result as
    ///necessary.
    class LIB_PUBLIC TermValarray : public std::valarray<double> {
        private:
            static unsigned __max_id;
#ifdef VERBOSE
            unsigned __id;
#endif
        public:
            TermValarray(size_t size = 0) :
                std::valarray<double>(size)
#ifdef VERBOSE
                , __id(__max_id++)
#endif
            {
#ifdef VERBOSE
                std::cout << "Constructed (size = " << size << ") Array "
                          << __id << std::endl;
#endif
            }

            TermValarray(double value, size_t size) :
                std::valarray<double>(value, size)
#ifdef VERBOSE
                , __id(__max_id++)
#endif
            {
#ifdef VERBOSE
                std::cout << "Constructed (size = " << size
                          << ", value =" << value
                          << ") Array " << __id
                          << std::endl;
#endif
            }

            TermValarray(const std::valarray<double> &orig) :
                std::valarray<double>(orig)
#ifdef VERBOSE
                , __id(__max_id++)
#endif
            {
#ifdef VERBOSE
                std::cout << "Constructed (copy) Array " << __id
                          << "from " << orig
                          << std::endl;
#endif
            }

            TermValarray(const double *data, size_t data_size) :
                std::valarray<double>(data, data_size)
#ifdef VERBOSE
                , __id(__max_id++)
#endif
            {
#ifdef VERBOSE
                std::cout << "Constructed (copy) Array " << __id
                          << "from double*: " << *this
                          << std::endl;
#endif
            }

            TermValarray &operator=(const std::valarray<double> &rhs)
            {
#ifdef VERBOSE
                std::cout << "Array " << __id
                          << " (size " << size()
                          << ") setting to " << rhs
                          << std::endl;
#endif
                if ( size() != rhs.size() ) resize(rhs.size());
                std::valarray<double>::operator=(rhs);
                return *this;
            }


#ifdef VERBOSE
            void resize(size_t size, double value = double())
            {
                std::cout << "Array " << __id
                          << " resizing to " << size << "x" << value
                          << std::endl;
                std::valarray<double>::resize(size, value);
            }
#endif

            TermValarray &operator*=(const std::valarray<double> &rhs)
            {
                std::valarray<double>::operator*=(rhs);
                return *this;
            }

            TermValarray &operator/=(const std::valarray<double> &rhs)
            {
                std::valarray<double>::operator/=(rhs);
                return *this;
            }

            TermValarray &operator+=(const std::valarray<double> &rhs)
            {
                std::valarray<double>::operator+=(rhs);
                return *this;
            }

            TermValarray &operator-=(const std::valarray<double> &rhs)
            {
                std::valarray<double>::operator-=(rhs);
                return *this;
            }
    };

#ifdef VERBOSE
    unsigned TermValarray::__max_id = 0;
#endif

    namespace TermCalculator {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        ///A class to hold the variables that may appear it PSF terms.
        template<typename ResultType>
        class Variables : public qi::symbols<char, ResultType> {
            public:
                Variables() {}
        };

        ///The grammar defining the calculator.
        template<typename Iterator, typename ResultType>
        class Grammar : public qi::grammar<Iterator,
                                           ResultType(),
                                           ascii::space_type> {
            private:
                ///The size of the resulting arrays.
                size_t __result_size;

                ///Alias of std::pow, needed for type determination.
                ResultType pow(const std::valarray<double> &base,
                               const std::valarray<double> &exponent)
                {return ResultType(std::pow(base, exponent));}

                ///Alias of std::cos, needed for type determination.
                ResultType cos(const std::valarray<double> &x)
                {return ResultType(std::cos(x));}

                ///Alias of std::sin, needed for type determination.
                ResultType sin(const std::valarray<double> &x)
                {return ResultType(std::sin(x));}

                ///Alias of std::tan, needed for type determination.
                ResultType tan(const std::valarray<double> &x)
                {return ResultType(std::tan(x));}

                ///Alias of std::acos, needed for type determination.
                ResultType acos(const std::valarray<double> &x)
                {return ResultType(std::acos(x));}

                ///Alias of std::asin, needed for type determination.
                ResultType asin(const std::valarray<double> &x)
                {return ResultType(std::asin(x));}

                ///Alias of std::atan, needed for type determination.
                ResultType atan(const std::valarray<double> &x)
                {return ResultType(std::atan(x));}

                ///Alias of std::atan2, needed for type determination.
                ResultType atan2(const std::valarray<double> &x,
                                 const std::valarray<double> &y)
                {return ResultType(std::atan2(x, y));}

                ///Alias of std::cosh, needed for type determination.
                ResultType cosh(const std::valarray<double> &x)
                {return ResultType(std::cosh(x));}

                ///Alias of std::sinh, needed for type determination.
                ResultType sinh(const std::valarray<double> &x)
                {return ResultType(std::sinh(x));}

                ///Alias of std::tanh, needed for type determination.
                ResultType tanh(const std::valarray<double> &x)
                {return ResultType(std::tanh(x));}

                ///Alias of std::exp, needed for type determination.
                ResultType exp(const std::valarray<double> &x)
                {return ResultType(std::exp(x));}

                ///Alias of std::log, needed for type determination.
                ResultType log(const std::valarray<double> &x)
                {return ResultType(std::log(x));}

                ///Alias of std::log10, needed for type determination.
                ResultType log10(const std::valarray<double> &x)
                {return ResultType(std::log10(x));}

                ///Alias of std::sqrt, needed for type determination.
                ResultType sqrt(const std::valarray<double> &x)
                {return ResultType(std::sqrt(x));}

                ///Alias of std::ceil, needed for type determination.
                ResultType ceil(const ResultType &x)
                {
                    ResultType result(x.size());
                    for( size_t i = 0; i < x.size(); ++i )
                        result[i]=std::ceil(x[i]);
                    return result;
                }

                ///Alias of std::floor, needed for type determination.
                ResultType floor(const ResultType &x)
                {
                    ResultType result(x.size());
                    for( size_t i = 0; i < x.size(); ++i )
                        result[i]=std::floor(x[i]);
                    return result;
                }

                ///Alias of std::fmod, needed for type determination.
                ResultType fmod(const ResultType &numer,
                                 const ResultType &denom)
                {
                    assert(numer.size() == denom.size());
                    ResultType result(numer.size());
                    for( size_t i = 0; i < numer.size(); ++i )
                        result[i]=std::fmod(numer[i], denom[i]);
                    return result;
                }

                ///Alias of std::abs, needed for type determination.
                ResultType abs(const std::valarray<double> &x)
                {return ResultType(std::abs(x));}


                qi::rule<Iterator, ResultType(), ascii::space_type>
                    ///Rule defining the final expression as a sum of terms.
                    __expression,

                    ///\brief Rule defining a term in the expression as a
                    ///product of factors.
                    __factor,

                    ///\brief Rule defining a value in a power as unary
                    ///function or number.
                    __value,

                    ///Rule defining exponentiation.
                    __exponentiation;


                ///\brief The user defined variables that may appear in the
                ///expressions.
                Variables<ResultType> __variables;

            public:

                ///\brief Construct a term calculator producing results of
                ///the given length.
                Grammar(size_t result_size=0);

                ///Define length of the final result vector.
                void set_result_size(size_t result_size)
                {__result_size=result_size;}

                ///\brief Define a new variable which may participate in the
                ///expressions evaluated.
                void add_variable(
                        ///The name to use for the variable.
                        const std::string &name,

                        ///The actual vector of values for the variable. The
                        ///size must match the currently defined result size.
                        const ResultType &value
                );

                ///Unmodifiable reference to the defined variables.
                const Variables<ResultType> &variables() const
                {return __variables;}
        };
    }

    namespace TermGenerator {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        template<typename ResultType>
        struct PolyOrderSpec {
            unsigned max_order;
            ResultType terms;
        };

        ///The grammar of the term generator language.
        template<typename Iterator, typename ResultType>
        class Grammar : public boost::spirit::qi::grammar<Iterator,
                                                          ResultType(),
                                                          ascii::space_type> {
            private:
                ///Rule defining a single term in a comma separated list.
                qi::rule<Iterator, std::string(), ascii::space_type>
                    __simple_term;

                qi::rule<Iterator, ResultType(), ascii::space_type>

                    ///\brief Rule defining a comma separated list of terms
                    ///in {...}.
                    __term_list,

                    ///Rule defining a set of terms.
                    __term_set,

                    ///\brief Rule defining a cross product of an arbitrary
                    ///number of term sets.
                    __cross_product,

                    ///\brief Rule defining the full generator as union of
                    ///cross products.
                    __expression;

                ///Add all terms of a giver order to the result.
                void add_poly_terms(
                    ///The order of the term to add (cumulative).
                    unsigned order,

                    ///The first term to use in the expansion.
                    const typename ResultType::const_iterator&
                    first_input_term,

                    ///One past the last term to use in the expansion.
                    const typename ResultType::const_iterator&
                    last_input_term,

                    ///The container to add the generated terms to.
                    ResultType &output_terms,

                    ///Prefix to prepend to each added term.
                    const std::string &prefix = ""
                );

                ///\brief Combine terms in all possible ways up to a maximum
                ///combined order.
                void order(
                    ///The maximum cumulative order of all terms involved.
                    unsigned max_order,

                    ///The terms to combine.
                    const ResultType &input_terms,

                    ///The target to fill with the combined terms.
                    ResultType &output_terms
                );

                ///\brief The cross product of two sets of terms.
                void cross_set(ResultType &lhs, const ResultType &rhs);

                ///Combine the terms of two sets.
                void union_set(ResultType &lhs, const ResultType &rhs)
                {lhs.insert(lhs.end(), rhs.begin(), rhs.end());}

            public:
                ///Create an instance.
                Grammar();

                ///Describe the currently supported language (EBNF format).
        };

        const std::string ebnf_definition =
            "(* items in angle brackets (< or >) are assumed to be obvious "
            "and thus not defined. *)\n"
            "termchar = <ascii character> - \",\" - \"}\" ;\n"
            "term = termchar , { termchar } ; "
            "(* mathematical expressions involving variables, floating point"
            " numbers and Pi. The complete list of mathematical functions "
            "from c++99's cmath library are supported. *)\n"
            "list = \"{\" , term , { \",\" , term } , \"}\"; "
            "(* simple listing of terms to include *)\n"
            "poly = \"O\" , <integer> , list; "
            "(* Expands to all polynomial terms of up to "
            "combined order <integer> of the entries in list *)\n"
            "set = list | poly ;"
            "cross = set , { \"*\" , set } ;"
            "(* expands to the cross product of all sets. *)\n"
            "expression = cross , { \"+\" , cross } (* merge the terms of "
            "all cross products together. *)\n";
    }

    namespace TermCalculator {

        template<typename Iterator, typename ResultType>
        Grammar<Iterator, ResultType>::Grammar(size_t result_size) :
            Grammar::base_type(__expression),
            __result_size(result_size)
        {
            using qi::double_;
            using qi::_val;
            using qi::_1;
            using boost::phoenix::bind;

            __exponentiation = (
                    __value[_val = _1]
                    >>
                    -('^' >> __value[_val = bind(&Grammar::pow,
                                                 this,
                                                 _val,
                                                 _1)])
            );

            __value = (
                    double_[bind(&ResultType::resize,
                                 &_val,
                                 __result_size,
                                 _1)]

                    | __variables[_val = _1]

                    | ( "cos("
                        >>
                        __expression[_val = bind(&Grammar::cos, this, _1)]
                        >>
                        ')' )

                    | ( "sin("
                        >>
                        __expression[_val = bind(&Grammar::sin, this, _1)]
                        >>
                        ')' )

                    | ( "tan("
                        >>
                        __expression[_val = bind(&Grammar::tan, this, _1)]
                        >>
                        ')' )

                    | ( "acos("
                        >>
                        __expression[_val = bind(&Grammar::acos, this, _1)]
                        >>
                        ')' )

                    | ( "asin("
                        >>
                        __expression[_val = bind(&Grammar::asin, this, _1)]
                        >>
                        ')' )

                    | ( "atan("
                        >>
                        __expression[_val = bind(&Grammar::atan, this, _1)]
                        >>
                        ')' )

                    | ( "atan2("
                        >>
                        __expression[_val = _1]
                        >>
                        ','
                        >>
                        __expression[
                            _val = bind(&Grammar::atan2, this, _val, _1)
                        ]
                        >>
                        ')' )

                    | ( "cosh("
                        >>
                        __expression[_val = bind(&Grammar::cosh, this, _1)]
                        >>
                        ')' )

                    | ( "sinh("
                        >>
                        __expression[_val = bind(&Grammar::sinh, this, _1)]
                        >>
                        ')' )

                    | ( "tanh("
                        >>
                        __expression[_val = bind(&Grammar::tanh, this, _1)]
                        >>
                        ')' )

                    | ( "exp("
                        >>
                        __expression[_val = bind(&Grammar::exp, this, _1)]
                        >>
                        ')' )

                    | ( "log("
                        >>
                        __expression[_val = bind(&Grammar::log, this, _1)]
                        >>
                        ')' )

                    | ( "log10("
                        >>
                        __expression[_val = bind(&Grammar::log10, this, _1)]
                        >>
                        ')' )

                    | ( "sqrt("
                        >>
                        __expression[_val = bind(&Grammar::sqrt, this, _1)]
                        >>
                        ')' )

                    | ( "ceil("
                        >>
                        __expression[_val = bind(&Grammar::ceil, this, _1)]
                        >>
                        ')' )

                    | ( "floor("
                        >>
                        __expression[_val = bind(&Grammar::floor, this, _1)]
                        >>
                        ')' )

                    | ( "fmod("
                        >>
                        __expression[_val = _1]
                        >>
                        ','
                        >>
                        __expression[_val = bind(&Grammar::fmod,
                                                 this,
                                                 _val,
                                                 _1)]
                        >>
                        ')' )

                    | ( "abs("
                        >>
                        __expression[_val = bind(&Grammar::abs, this, _1)]
                        >>
                        ')' )

                    | ( '(' >> __expression[_val = _1] >> ')' )
            );

            __factor = (
                    __exponentiation[_val = _1]
                    >>
                    *(
                        ( '*' >> __exponentiation[_val *= _1] )
                        | ( '/' >> __exponentiation[_val /= _1] )
                    )
            );

            __expression = (
                    __factor[_val = _1]
                    >>
                    *(
                        ( '+' >> __factor[_val += _1] )
                        |
                        ( '-' >> __factor[_val -= _1] )
                    )
            );
        }

        template<typename Iterator, typename ResultType>
        void Grammar<Iterator, ResultType>::add_variable(
            const std::string &name,
            const ResultType &value
        )
        {
#ifdef VERBOSE
            std::cout << "Adding variable " << name << ": " << value
                      << std::endl;
#endif
            assert(value.size() == __result_size);
            __variables.add(name, value);
#ifdef VERBOSE
            std::cout << "Added variable " << name << ": " << value
                      << std::endl;
#endif
        }
    }

    namespace TermGenerator {
        namespace phoenix = boost::phoenix;
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        template<typename Iterator, typename ResultType>
        void Grammar<Iterator, ResultType>::add_poly_terms(
                unsigned order,
                const typename ResultType::const_iterator& first_input_term,
                const typename ResultType::const_iterator& last_input_term,
                ResultType &output_terms,
                const std::string &prefix)
        {
            typename ResultType::const_iterator
                next_input_term=first_input_term;
            ++next_input_term;
            if ( next_input_term == last_input_term ) {
                if ( order == 0 && prefix == "" )
                    output_terms.push_back("1");
                else {
                    std::ostringstream term_to_add;
                    if ( prefix != "" ) {
                        term_to_add << prefix;
                        if ( order != 0 ) term_to_add << " * ";
                    }
                    if ( order != 0 ) {
                        term_to_add << "(" << *first_input_term << ")";
                        if ( order != 1 )
                            term_to_add << "^" << order;
                    }
                    output_terms.push_back(term_to_add.str());
                }
            } else {
                for ( unsigned remaining_order = 0;
                        remaining_order <= order;
                        ++remaining_order ) {
                    std::ostringstream new_prefix;
                    new_prefix << prefix;
                    if ( remaining_order < order ) {
                        if ( prefix != "" ) new_prefix << " * ";
                        new_prefix << "(" + *first_input_term + ")";
                        if ( order - remaining_order > 1 )
                            new_prefix << "^" << order - remaining_order;
                    }
                    add_poly_terms(remaining_order,
                                   next_input_term,
                                   last_input_term,
                                   output_terms,
                                   new_prefix.str());
                }
            }
        }

        template<typename Iterator, typename ResultType>
        void Grammar<Iterator, ResultType>::order(unsigned max_order,
                                                  const ResultType &input_terms,
                                                  ResultType &output_terms)
        {
            output_terms.clear();
            for ( unsigned order = 0; order <= max_order; ++order )
                add_poly_terms(order,
                               input_terms.begin(),
                               input_terms.end(),
                               output_terms);
        }

        template<typename Iterator, typename ResultType>
        void Grammar<Iterator, ResultType>::cross_set(ResultType &lhs,
                                                      const ResultType &rhs)
        {
            ResultType result;
            for( typename ResultType::const_iterator
                    lhs_iter = lhs.begin();
                    lhs_iter != lhs.end();
                    ++lhs_iter )
                for ( typename ResultType::const_iterator
                        rhs_iter = rhs.begin();
                        rhs_iter != rhs.end();
                        ++rhs_iter )
                    result.push_back(
                            "(" + *lhs_iter + ") * (" + *rhs_iter + ")"
                    );
            lhs=result;
        }

        template<typename Iterator, typename ResultType>
        Grammar<Iterator, ResultType>::Grammar() :
            Grammar::base_type(__expression)
        {
            using phoenix::push_back;
            using phoenix::ref;
            using phoenix::at_c;
            using phoenix::bind;

            using qi::char_;
            using qi::uint_;
            using qi::_val;
            using qi::_1;
            using qi::_2;

            __simple_term = +( char_ - ',' - '}');

            __term_list = (
                    '{'
                    >> ( __simple_term[push_back(_val, _1)] % ',' )
                    >> '}'
            );

            __term_set = (
                    __term_list[_val = _1]
                    | ( 'O' >> uint_ >> __term_list )[
                          bind(&Grammar::order, this, _1, _2, _val)
                      ]
            );

            __cross_product = (
                    __term_set[_val = _1]
                    >> *( '*'
                          >>
                          __term_set[
                              bind(&Grammar::cross_set, this, _val, _1)
                          ] )
            );

            __expression = (
                    __cross_product[
                        bind(&Grammar::union_set, this, _val, _1)
                    ] % '+'
            );

        }
    }

    ///Fill result with the value of each term for each source.
    template<class TermIterator, class VarIterator, class ResultType>
        void evaluate_terms(
            ///Iterator to the first term to evaluate. Each term is a string
            ///expression involving variables defined by start_var/end_var.
            TermIterator start_term,

            ///One past the last term to evaluate.
            TermIterator end_term,

            ///See same name argument to evaluate_term_expression.
            VarIterator start_var,

            ///See same name argument to evaluate_term_expression.
            VarIterator end_var,

            ///See same name argument to evaluate_term_expression.
            ResultType &result
        )
    {
        TermCalculator::Grammar< std::string::const_iterator, TermValarray >
            calculator_grammar(start_var->second.size());

        for(; start_var != end_var; ++start_var) {
            calculator_grammar.add_variable(start_var->first,
                                            start_var->second);
        }

        unsigned result_index = 0;
        for(; start_term != end_term; ++start_term) {
            using boost::spirit::ascii::space;
            using boost::spirit::qi::phrase_parse;

            std::string::const_iterator parse_start = start_term->begin(),
                                        parse_end = start_term->end();
            bool parser_status = phrase_parse(parse_start,
                                              parse_end,
                                              calculator_grammar,
                                              space,
                                              result[result_index++]);
            if( !parser_status || parse_start != parse_end ) {
                std::ostringstream msg;
                msg << "Parsing " << *start_term
                    << " failed at " << std::string(parse_start, parse_end);
                throw Error::ParsingError(msg.str());
            }
        }
    }

    ///Fill result with all terms the PSF depends on per term_expression
    template<class VarIterator, class ResultType>
        void evaluate_term_expression(
            ///The expression defining the terms to include in PSF fitting (see
            ///Grammar for syntax).
            const std::string &term_expression,

            ///Iterator to the first varuable that participates in the
            ///expression. Entries should be pairs of \<variable name\>,
            ///\<array of values for each source\>
            VarIterator start_var,

            ///Iterator to one past the last variable that participates in the
            ///expression.
            VarIterator end_var,

            ///The location to fill with the PSF map terms. Resized as necessary
            ///to fit the result.
            ResultType &result
        )
    {
        assert(term_expression != "");
        using boost::spirit::ascii::space;
        using boost::spirit::qi::phrase_parse;
        std::string::const_iterator parse_start=term_expression.begin(),
                                    parse_end=term_expression.end();

        typedef std::list< std::string > GeneratorResultType;
        GeneratorResultType terms;
        TermGenerator::Grammar< std::string::const_iterator,
                                GeneratorResultType > generator_grammar;

        bool parser_status = phrase_parse(parse_start,
                                          parse_end,
                                          generator_grammar,
                                          space,
                                          terms);

        if ( !parser_status || parse_start != parse_end ) {
            std::ostringstream message;
            message << "PSF terms parsing failed at "
                    << std::string(parse_start, parse_end);
            throw Error::CommandLine(message.str());
        }

        result.resize(terms.size());
        evaluate_terms(terms.begin(),
                       terms.end(),
                       start_var,
                       end_var,
                       result);
    }

} //End PSF namespace.

#endif
