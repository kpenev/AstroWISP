#include "OutputArray.h"

namespace IO {
    template<>
        void OutputArray<double>::parse(const boost::any &value)
        {
            if(value.type() == typeid(double)) {
                __size = 1;
                __allocated_data = new double[1];
                __allocated_data[0] = boost::any_cast<const double&>(value);
                __data = __allocated_data;
                return;
            }
            if(value.type() == typeid(Core::RealList)) {
                const Core::RealList &list = 
                    TranslateToAny< Core::RealList >().get_value(value);
                __allocated_data = new double[list.size()];
                unsigned index = 0;
                for(
                    Core::RealList::const_iterator i = list.begin();
                    i != list.end();
                    ++i
                )
                    __allocated_data[index++] = *i;
                __size = list.size();
                __data = __allocated_data;
                return;
            }
            try {
                const std::vector<double>& vector =
                    TranslateToAny< std::vector<double> >().get_value(value);
                __size = vector.size();
                __data = &(vector[0]);
                return;
            } catch(boost::bad_any_cast) {}
            try {
                const std::valarray<double>& valarray =
                    TranslateToAny< std::valarray<double> >().get_value(
                        value
                    );
                __size = valarray.size();
                __data = &(valarray[0]);
                return;
            } catch(boost::bad_any_cast) {}
            try {
                const Eigen::VectorXd& 
                    vector = TranslateToAny< Eigen::VectorXd >().get_value(
                        value
                    );
                __size = vector.size();
                __data = vector.data();
                return;
            } catch(boost::bad_any_cast) {}
            const Eigen::ArrayXd& 
                array = TranslateToAny< Eigen::ArrayXd >().get_value(value);
            __size = array.size();
            __data = array.data();
        }

}
