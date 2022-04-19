#include "PSF.h"
#include "../Core/NaN.h"
#include <limits>
#include <algorithm>
#include <sstream>
#include <iostream>

namespace PSF {

    double PSF::line_circle_intersections(
            double x1, double x2, double y, double r, bool left_on_circle,
            bool right_on_circle) const
    {
        if(y>r) return 0;
        if(left_on_circle && right_on_circle) return Core::NaN;
        double r_2=r*r, y_2=y*y, left_dist_2=x1*x1+y_2,
               right_dist_2=x2*x2+y_2;
        if(left_on_circle) return (right_dist_2>r_2 ? 0 : Core::NaN);
        if(right_on_circle) return (left_dist_2>r_2 ? 0 : Core::NaN);
        if(left_dist_2<r_2 && right_dist_2<r_2) return Core::NaN;
        if(left_dist_2>r_2 && right_dist_2>r_2) {
            if(x1*x2>0) return 0;
            else return -std::sqrt(r_2-y_2);
        } else return std::sqrt(r_2-y_2);
    }

    double PSF::integrate_overlap_bottom_out(
        double x1,
        double y1,
        double x2,
        double y2,
        double rc,
        std::valarray<bool> on_circle,
        std::valarray<double> &intersections
    ) const
    {
        if(y1>=0 || intersections[2]==0) return 0.0;//entirely outside circ
        else if(on_circle[2] || on_circle[3] ||
                !std::isnan(intersections[2])) {//2 sides entirely outside
            if(x2>0) return (on_circle[3] ? 0.0 :
                    integrate_wedge(x1, y2, rc, false, true));
            else return (on_circle[2] ? 0.0 :
                    integrate_wedge(x2, y2, rc, true, true));
        } else { //1 side entirely out 1 side entirely in
            double y_split=std::numeric_limits<double>::infinity();
            if(!std::isnan(intersections[1])) y_split=intersections[1];
            if(!std::isnan(intersections[3]))
                    y_split=std::min(y_split, intersections[3]);
            y_split*=-1;
            return (integrate_rectangle((x1+x2)/2, (y2+y_split)/2,
                                        x2-x1, y2-y_split)+
                    (x2>0 ? integrate_wedge(x1, y_split, rc, false, true) :
                            integrate_wedge(x2, y_split, rc, true, true)));
        }
    }

    double PSF::integrate_overlap_bottom_in(
        double x1,
        double y1,
        double x2,
        double y2,
        double rc,
        std::valarray<bool> on_circle,
        std::valarray<double> &intersections
    ) const
    {
        double result=0;
        if(y1<0 || (std::isnan(intersections[1]) && x2>0) ||
           (std::isnan(intersections[3]) && x1<0)) {//entirely inside circle
            double result = integrate_rectangle((x1 + x2) / 2,
                                                (y1 + y2) / 2,
                                                x2 - x1,
                                                y2 - y1);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Rectangle integral over (" << x1 << ", " << y1
                    << "), (" << x2 << ", " << y2 << ") returned NaN";
                throw Error::Runtime(msg.str());
            }
            return result;
        }
        else if (!on_circle[0] && !on_circle[1]) {
            double y_split = y1;
            if(!std::isnan(intersections[1])) {
                if(x2>0) y_split=intersections[1];
                else {
                    y2=intersections[1];
                    on_circle[2]=true;
                }
            }
            if(!std::isnan(intersections[3])) {
                if(x2>0) {
                    y2=intersections[3];
                    on_circle[3]=true;
                } else y_split=intersections[3];
            }
            result=integrate_rectangle((x1+x2)/2, (y1+y_split)/2, x2-x1,
                                       y_split-y1);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Rectangle integral over (" << x1 << ", " << y1
                    << "), (" << x2 << ", " << y_split << ") returned NaN";
                throw Error::Runtime(msg.str());
            }
            y1=y_split;
            on_circle[(x2>0 ? 1 : 0)]=true;
        }
        if(x2>0 && !on_circle[3] && intersections[2]!=0) {
            result+=integrate_rectangle((x1+intersections[2])/2, (y1+y2)/2,
                                        intersections[2]-x1, y2-y1);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Rectangle integral over (" << x1 << ", " << y1
                    <<" ), (" << intersections[2] << ", " << y2
                    << ") returned NaN";
                throw Error::Runtime(msg.str());
            }
            on_circle[3]=true;
            x1=intersections[2];
            intersections[2]=0;
        } else if(x1<0 && !on_circle[2] && intersections[2]!=0) {
            result+=integrate_rectangle((x2-intersections[2])/2, (y1+y2)/2,
                                        x2+intersections[2], y2-y1);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Rectangle integral over (" << x2 << ", " << y1
                    <<" ), (" << intersections[2] << ", " << y2
                    << ") returned NaN";
                throw Error::Runtime(msg.str());
            }
            on_circle[2]=true;
            x2=-intersections[2];
            intersections[2]=0;
        }
        if(x2>0 && !on_circle[0]) {
            result+=integrate_wedge(x1, y1, rc, false, false);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Wedge integral over (" << x1 << ", " << y1 << ", "
                    << rc << ", right, top) returned NaN";
                throw Error::Runtime(msg.str());
            }
        } else if(!on_circle[1]) {
            result+=integrate_wedge(x2, y1, rc, true, false);
            if(std::isnan(result)) {
                std::ostringstream msg;
                msg << "Wedge integral over (" << x1 << ", " << y1 << ", "
                    << rc << ", left, top) returned NaN";
                throw Error::Runtime(msg.str());
            }
        }
        return result;
    }

    double PSF::integrate_overlap(
            double x1, double y1, double x2, double y2, double rc,
            const std::valarray<bool> &on_circle) const
    {
        if(x1*x2<0) {
            std::valarray<bool> on_circle_copy(on_circle);
            on_circle_copy[1]=on_circle_copy[2]=false;
            double result=integrate_overlap(x1,
                                            y1,
                                            0,
                                            y2,
                                            rc,
                                            on_circle_copy);
            on_circle_copy[1]=on_circle[1];
            on_circle_copy[2]=on_circle[2];
            on_circle_copy[0]=on_circle_copy[3]=false;
            return result+integrate_overlap(0,
                                            y1,
                                            x2,
                                            y2,
                                            rc,
                                            on_circle_copy);
        }
        if(y1*y2<0) {
            std::valarray<bool> on_circle_copy(on_circle);
            on_circle_copy[2]=on_circle_copy[3]=false;
            double result=integrate_overlap(x1,
                                            y1,
                                            x2,
                                            0,
                                            rc, on_circle_copy);
            on_circle_copy[2]=on_circle[2];
            on_circle_copy[3]=on_circle[3];
            on_circle_copy[0]=on_circle_copy[1]=false;
            return result+integrate_overlap(x1,
                                            0,
                                            x2,
                                            y2,
                                            rc,
                                            on_circle_copy);
        }

        //Now entire rectangle is restricted to a single quadrant
        std::valarray<double> intersections(4);
        intersections[0]=line_circle_intersections(x1, x2, std::abs(y1), rc,
                on_circle[0], on_circle[1]);
        intersections[1]=line_circle_intersections(y1, y2, std::abs(x2), rc,
                        on_circle[1], on_circle[2]);
        intersections[2]=line_circle_intersections(x1, x2, std::abs(y2), rc,
                        on_circle[3], on_circle[2]);
        intersections[3]=line_circle_intersections(y1, y2, std::abs(x1), rc,
                        on_circle[0], on_circle[3]);
        if(intersections[0]==0)
            return integrate_overlap_bottom_out(x1,
                                                y1,
                                                x2,
                                                y2,
                                                rc,
                                                on_circle,
                                                intersections);
        else if(std::isnan(intersections[0]))
            return integrate_overlap_bottom_in(x1, y1, x2, y2, rc, on_circle,
                    intersections);
        else {
            std::valarray<bool> on_circle_copy(on_circle);
            if(y2>0) {
                if(x2>0) {
                    on_circle_copy[1]=true;
                    x2=intersections[0];
                    intersections[0]=Core::NaN;
                } else {
                    on_circle_copy[0]=true;
                    x1=-intersections[0];
                    intersections[0]=Core::NaN;
                }
                return integrate_overlap_bottom_in(x1, y1, x2, y2, rc,
                        on_circle_copy, intersections);
            } else {
                double result=0;
                if(x2>0) {
                    result=integrate_rectangle((x1+intersections[0])/2,
                            (y1+y2)/2, intersections[0]-x1, y2-y1);
                    x1=intersections[0];
                    intersections[0]=0;
                    on_circle_copy[0]=true;
                } else {
                    result=integrate_rectangle((x2-intersections[0])/2,
                            (y1+y2)/2, x2+intersections[0], y2-y1);
                    x2=-intersections[0];
                    intersections[0]=0;
                    on_circle_copy[1]=true;
                }
                return result+integrate_overlap_bottom_out(x1,
                                                           y1,
                                                           x2,
                                                           y2,
                                                           rc,
                                                           on_circle_copy,
                                                           intersections);
            }
        }
    }

    double PSF::integrate(
        double center_x,
        double center_y,
        double dx,
        double dy,
        double circle_radius
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
        ,
        bool reset_piece_id,
        bool skip_piece
#endif
#endif
    ) const
    {
        if(circle_radius) {
            double dxh=dx/2, dyh=dy/2;
            return integrate_overlap(center_x-dxh, center_y-dyh, center_x+dxh,
                        center_y+dyh, circle_radius);
        } else return integrate_rectangle(center_x, center_y, dx, dy);
    }

} //End PSF namespace.
