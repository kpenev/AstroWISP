#include "WedgeIntegral.h"
#include "../Core/NaN.h"

namespace PSF {

    void WedgeIntegral::fill_p_even(Core::vector_size_type max_n)
    {
        if(__deltay_pow.size() <= max_n + 1)
            fill_powers(__deltay_pow, max_n + 1);
        for(Core::vector_size_type n = __q[0].size(); n <= max_n; n += 2) {
            __q[0].push_back(__deltay_pow[n + 1] * 2.0 / (n + 1));
            __q[0].push_back(0);
        }
    }

    void WedgeIntegral::fill_p_odd(Core::vector_size_type max_n)
    {
        if(__deltay_pow.size() < max_n) fill_powers(__deltay_pow, max_n - 1);

        double p2 = __q[1][__q[1].size() - 2],
               p1 = __q[1][__q[1].size() - 1];
        for(Core::vector_size_type n = __q[1].size(); n <= max_n; ++n) {
            double new_podd = (__r2_m_y02 * (n - 1) * p2
                               -
                               __y0 * (2 * n + 1) * p1
                               +
                               __deltay_pow[n - 1] * (
                                   n % 2
                                   ? __p_odd_free_term_mult_odd
                                   : __p_odd_free_term_mult_even
                               )
            ) / (2.0 + n);
            __q[1].push_back(new_podd);
            p2 = p1;
            p1 = new_podd;
        }
    }

    void WedgeIntegral::fill_q_diagonal(Core::vector_size_type m,
                                        Core::vector_size_type n)
    {
        Core::vector_size_type m_fill = m, n_fill = n;
        while(std::isnan(__q[m_fill][n_fill])) {
            m_fill -= 2;
            n_fill += 2;
        }
        for(;m_fill < m; m_fill += 2, n_fill -= 2)
            __q[m_fill + 2][n_fill - 2] = (
                __r2_m_y02 * __q[m_fill][n_fill - 2]
                -
                __twice_y0 * __q[m_fill][n_fill - 1]
                -
                __q[m_fill][n_fill]
            );
    }

    void WedgeIntegral::calculate_q(Core::vector_size_type m,
                                    Core::vector_size_type n)
    {
        if(m % 2) {
            if(__q[1].size() < m + n + 2) fill_p_odd(m + n - 1);
        } else if(__q[0].size() <= m + n) fill_p_even(m + n);
        Core::vector_size_type m_fill = m;
        while(std::isnan(__q[m_fill][n])) m_fill -= 2;
        for(;m_fill < m; m_fill += 2) {
            if(std::isnan(__q[m_fill][n + 1]))
                fill_q_diagonal(m_fill, n + 1);
            fill_q_diagonal(m_fill + 2, n);
        }
    }

#ifdef DEBUG
    void WedgeIntegral::output_q(std::ostream &os)
    {
        os.precision(10);
        os.setf(std::ios_base::scientific);
        for(int n = std::max(__q[0].size(), __q[1].size()) - 1; n >= 0; --n){
            for(
                int m = 0;
                !std::isnan(__q[m][n]) || !std::isnan(__q[m+1][n]);
                ++m
            )
                os << std::setw(20) << __q[m][n];
            os << std::endl;
        }
    }
#endif

    void WedgeIntegral::expand_storage(Core::vector_size_type min_m_size,
                                       Core::vector_size_type min_n_size)
    {
        Core::vector_size_type m_size = 2 * __values.size(),
                               n_size = 2 * __values[0].size();
        while(m_size < min_m_size) m_size *= 2;
        while(n_size < min_n_size) n_size *= 2;
        __values.resize(m_size);
        __q.resize(m_size);
        for(Core::vector_size_type m = 0; m < m_size; ++m) {
            __values[m].resize(n_size, Core::NaN);
            if(m > 1) __q[m].resize(n_size, Core::NaN);
        }
    }

    WedgeIntegral::WedgeIntegral(double r, double x0, double y0) :
        __x0(x0),
        __y0(y0),
        __y02(y0 * y0),
        __twice_y0(2.0 * y0),
        __r2(r * r),
        __r2_m_y02(__r2 - __y02),
        __values(__initial_storage),
        __q(__initial_storage)
    {
#ifdef DEBUG
        assert(x0 > 0 && y0 > 0 && r > 0);
#endif
        initialize_powers(__x0_pow, x0, 3, __initial_storage);

        double x02_p_y02 = __x0_pow[2] + __y02,
               delta_factor;
        if(__r2 <= x02_p_y02) delta_factor = 0;
        else delta_factor = std::sqrt((__r2 - x02_p_y02) / x02_p_y02);
        __deltax = y0 * delta_factor;
        __deltay = x0 * delta_factor;
        __x_min = std::max(x0 - __deltax, 0.0);
        __x_max = x0 + __deltax;
        __y_min = std::max(y0 - __deltay, 0.0);
        __y_max = y0 + __deltay;

        initialize_powers(__deltax_pow, __deltax, 3, __initial_storage);
        initialize_powers(__deltay_pow, __deltay, 1, __initial_storage);

        __p_odd_free_term_mult_even = -2.0 * (__x0_pow[3]
                                              +
                                              3.0 * __x0 * __deltax_pow[2]);
        __p_odd_free_term_mult_odd = 2.0 * (3.0 * __x0_pow[2] * __deltax
                                            +
                                            __deltax_pow[3]);

        for(Core::vector_size_type i = 0; i < __initial_storage; ++i) {
            __values[i].resize(__initial_storage, Core::NaN);
            if(i > 1) __q[i].resize(__initial_storage, Core::NaN);
        }
        __q[0].reserve(__initial_storage);
        __q[1].reserve(__initial_storage);


        double x_max3 = std::pow(__x_max, 3),
               x_min3 = std::pow(__x_min, 3),
               p_odd_common = 0.5 * (
                   (__y_max * __x_min - __y_min * __x_max)
                   +
                   __r2 * (
                       (__y_min > 0 ? std::atan(__x_max / __y_min) : M_PI/2)
                       -
                       std::atan(__x_min / __y_max)
                   )
               );
        __q[1].push_back(p_odd_common);
        __q[1].push_back((x_max3 - x_min3) / 3 - y0 * p_odd_common);
        fill_p_even(1);

        __values[0][0] = (__q[1][0]
                          -
                          __q[0][0] * x0
                          +
                          2.0 * __deltax * __deltay);
    }

    double WedgeIntegral::operator()(unsigned m, unsigned n)
    {
        unsigned mp1 = m + 1,
                 np1 = n + 1;
        if(mp1 >= __values.size() || mp1 + np1>=__values[0].size())
            expand_storage(mp1, mp1 + np1);
        if(std::isnan(__values[m][n])) {
            if(__deltax_pow.size() <= m + 1) fill_powers(__deltax_pow, mp1);
            if(__deltay_pow.size() <= n + 1) fill_powers(__deltay_pow, np1);
            if(__x0_pow.size() <= m + 1) fill_powers(__x0_pow, mp1);
            if(n % 2) __values[m][n] = 0;
            else {
                __values[m][n] = (2.0 * __deltax_pow[mp1] * __deltay_pow[np1]
                                  /
                                  np1);
                if(m % 2) __values[m][n] = -__values[m][n];
            }
            unsigned mp1_choose_i = 1;
            int sign=(m % 2 ? 1 : -1);
            for(unsigned i = 0; i <= mp1; ++i) {
                if(i == 0 && __q[0].size() <= n) fill_p_even(n);
                else if(i == 1 && __q[1].size() <= n) fill_p_odd(n);
                else if(std::isnan(__q[i][n])) calculate_q(i, n);
                __values[m][n] += (mp1_choose_i
                                   *
                                   __q[i][n]
                                   *
                                   __x0_pow[mp1 - i]
                                   *
                                   sign);
                sign = -sign;
                mp1_choose_i *= (mp1 - i);
                mp1_choose_i /= (i + 1);
            }
            __values[m][n] /= mp1;
        }
        return __values[m][n];
    }

} //End PSF namespace.
