#include "CirclePieceIntegral.h"
#include "../Core/NaN.h"

namespace PSF {

    void CirclePieceIntegral::fill_R(Core::vector_size_type max_n,
            std::vector<double> &y_pow, std::vector<double> &Rn)
    {
        if(y_pow.size()<=max_n) fill_powers(y_pow, max_n);
        double r0=Rn[0];
        for(Core::vector_size_type n=Rn.size(); n<=max_n; ++n)
            Rn.push_back(r0*y_pow[n]);
    }

    void CirclePieceIntegral::fill_p_even(Core::vector_size_type max_n)
    {
        if(__ymin_pow.size()<=max_n+1) fill_powers(__ymin_pow, max_n+1);
        if(__ymax_pow.size()<=max_n+1) fill_powers(__ymax_pow, max_n+1);
        for(Core::vector_size_type n=__q[0].size(); n<=max_n; ++n)
            __q[0].push_back((__ymax_pow[n+1]-__ymin_pow[n+1])/(n+1));
    }

    void CirclePieceIntegral::fill_p_odd(Core::vector_size_type max_n)
    {
        if(__R_ymin.size()<max_n) fill_R(max_n-1, __ymin_pow, __R_ymin);
        if(__R_ymax.size()<max_n) fill_R(max_n-1, __ymax_pow, __R_ymax);

        double p2=__q[1][__q[1].size()-2], p1=__q[1][__q[1].size()-1];
        for(Core::vector_size_type n=__q[1].size(); n<=max_n; ++n) {
            double new_podd=(
                    __R_ymin[n-1]-__R_ymax[n-1]
                    +
                    (2.0*n+1.0)*__yc*p1
                    +
                    __r2_m_yc2*(n-1)*p2)/(2.0+n);
            __q[1].push_back(new_podd);
            p2=p1;
            p1=new_podd;
        }
    }

    void CirclePieceIntegral::fill_q_diagonal(Core::vector_size_type m,
                                              Core::vector_size_type n)
    {
        Core::vector_size_type m_fill=m, n_fill=n;
        while(std::isnan(__q[m_fill][n_fill])) {m_fill-=2; n_fill+=2;}
        double twice_yc=2.0*__yc;
        for(;m_fill<m; m_fill+=2, n_fill-=2)
            __q[m_fill+2][n_fill-2]=
                __r2_m_yc2*__q[m_fill][n_fill-2]
                +
                twice_yc*__q[m_fill][n_fill-1]
                -
                __q[m_fill][n_fill];
    }

    void CirclePieceIntegral::calculate_q(Core::vector_size_type m,
                                          Core::vector_size_type n)
    {
        if(m%2) {
            if(__q[1].size()<m+n+2) fill_p_odd(m+n-1);
        } else if(__q[0].size()<=m+n) fill_p_even(m+n);
        Core::vector_size_type m_fill=m;
        while(std::isnan(__q[m_fill][n])) m_fill-=2;
        for(;m_fill<m; m_fill+=2) {
            if(std::isnan(__q[m_fill][n+1])) fill_q_diagonal(m_fill, n+1);
            fill_q_diagonal(m_fill+2, n);
        }
    }

#ifdef DEBUG
    void CirclePieceIntegral::assert_parameter_consistency()
    {
        assert(__ymin_pow[1]<=__ymax_pow[1]);
        assert(__xmin_pow[1]>=__xc_pow[1]);
        assert(__ymin_pow[1]>=0);
    }
#endif

    CirclePieceIntegral::CirclePieceIntegral(
        double xmin,
        double ymin,
        double ymax,
        double xc,
        double yc,
        double r,
        Core::vector_size_type initial_storage
    ) :
        __yc(yc),
        __r2(r * r),
        __yc2(yc * yc),
        __r2_m_yc2(__r2 - __yc2),
        __ymin_m_yc(ymin - yc),
        __ymax_m_yc(ymax - yc),
        __root_ymin(
            std::max(0.0, std::sqrt(__r2 - std::pow(__ymin_m_yc, 2)))
        ),
        __root_ymax(
            std::max(0.0, std::sqrt(__r2 - std::pow(__ymax_m_yc, 2)))
        ),
        __values(++initial_storage),
        __q(initial_storage + 1)
    {
        initialize_powers(__ymin_pow, ymin, 2, initial_storage);
        initialize_powers(__ymax_pow, ymax, 2, initial_storage);
        initialize_powers(__xc_pow, xc, 2, initial_storage);
        initialize_powers(__xmin_pow, xmin, 2, initial_storage);

#ifdef DEBUG
        assert_parameter_consistency();
#endif 

        __R_ymin.reserve(initial_storage);
        __R_ymin.push_back(std::pow(__root_ymin, 3));
        __R_ymax.reserve(initial_storage);
        __R_ymax.push_back(std::pow(__root_ymax, 3));

        Core::vector_size_type q_size=2*initial_storage+2;
        for(Core::vector_size_type i = 0; i < initial_storage; ++i) {
            __values[i].resize(initial_storage, Core::NaN);
            if(i) __q[i + 1].resize(q_size, Core::NaN);
        }
        __q[0].reserve(q_size);
        __q[1].reserve(q_size);

        double pie_area = __r2 / 2.0 * (
            std::atan(__ymax_m_yc / __root_ymax)
            -
            std::atan(__ymin_m_yc / __root_ymin)
        );
        __q[1].push_back(
            (__root_ymax * __ymax_m_yc - __root_ymin*__ymin_m_yc) / 2.0
            +
            pie_area
        );
        double twice_r2_p_yc2 = 2.0 * __r2 + __yc2;
        __q[1].push_back(
                (
                    __root_ymin * (twice_r2_p_yc2
                                   +
                                   yc*ymin
                                   -
                                   2.0 * __ymin_pow[2])
                    -
                    __root_ymax * (twice_r2_p_yc2
                                   +
                                   yc*ymax
                                   -
                                   2.0*__ymax_pow[2])
                ) / 6.0
                +
                __yc * pie_area
        );

        fill_p_even(1);

        __values[0][0] = __q[1][0] + __q[0][0] * (xc - xmin);
    }

    double CirclePieceIntegral::operator()(unsigned m, unsigned n)
    {
        unsigned mp1=m+1;
        if(m>=__values.size() || n>=__values[m].size())
            throw Error::Runtime(
                    "Resizing CirclePieceIntegrals not implemented");
        if(std::isnan(__values[m][n])) {
            fill_p_even(n);
            if(__xc_pow.size()<=m+1) fill_powers(__xc_pow, mp1);
            if(__xmin_pow.size()<=m+1) fill_powers(__xmin_pow, mp1);
            __values[m][n]=-__q[0][n]*__xmin_pow[mp1];
            unsigned mp1_choose_i=1;
            if(__q[1].size()<=n) fill_p_odd(n);
            for(unsigned i=0; i<=mp1; ++i) {
                if(i>1 && std::isnan(__q[i][n])) calculate_q(i, n);
                __values[m][n]+=mp1_choose_i*__q[i][n]*__xc_pow[mp1-i];
                mp1_choose_i*=(mp1-i);
                mp1_choose_i/=(i+1);
            }
            __values[m][n]/=mp1;
        }
        return __values[m][n];
    }

}
