#include "IntegralUtil.h"

namespace PSF {

    void fill_powers(std::vector<double> &powers,
                     Core::vector_size_type max_pow)
    {
        double x=powers[1], xn=powers.back();
        for(Core::vector_size_type n=powers.size(); n<=max_pow; n++) {
            xn*=x;
            powers.push_back(xn);
        }
    }

    void initialize_powers(std::vector<double> &powers,
                           double x,
                           Core::vector_size_type max_pow,
                           Core::vector_size_type initial_storage)
    {
        powers.reserve(std::max(initial_storage, max_pow));
        powers.push_back(1);
        powers.push_back(x);
        if(max_pow>1) fill_powers(powers, max_pow);
    }
}
