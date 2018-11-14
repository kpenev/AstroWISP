#include "OutputArray.h"

namespace IO {
    template<>
        LIB_PUBLIC void OutputArray<double>::parse(const boost::any &value)
        {
            if(value.type() == typeid(double)) {
                __size = 1;
                __allocated_data = new double[1];
                __allocated_data[0] = boost::any_cast<const double&>(value);
                __data = __allocated_data;
                return;
            }
            if(value.type() == typeid(Core::RealList)) {
                if(!try_container_type< Core::RealList >(value))
                    throw boost::bad_any_cast();
                return;
            }
            if(
                !(
                    try_container_type< std::vector<double> >(value)
                    ||
                    try_array_type< std::valarray<double> >(value)
                    ||
                    try_array_type< Eigen::VectorXd >(value)
                    ||
                    try_array_type< Eigen::ArrayXd >(value)
                )
            )
                throw boost::bad_any_cast();
        }

}
