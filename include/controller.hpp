#ifndef CONTROLLER_DEMO_H_
#define CONTROLLER_DEMO_H_

#include <iostream>
#include <vector>
#include "config.h"

class Controller {
public:
    Controller(double dt, double theta0) : dt_(dt), theta0_(theta0), theta_(theta0) {}
    ~Controller(){}

    double run(double yr, double x1, double x2, double x3, 
               double F, double dyr, double ddyr, double dddyr) {
        // update
        // updateYrDeri(yr)
        dyr_ = dyr;
        ddyr_ = ddyr;
        dddyr_ = dddyr;

        // calculation
        const double z1 = calcZ1(x1, yr);
        const double alpha1 = clacAlpha1(z1, dyr_);
        // const double x2 = calcX2(x1)
        const double z2 = calcZ2(x2, alpha1);
        const double alpha2 = calcAlpha2(z1, z2, x2, F, dyr_, ddyr_);
        const double z3 = calcZ3(x3, alpha2);

        updateTheta(z3, x3, theta_);

        const double u = calcU(z2, z3, x2, theta_, x3, v1_, v2_, v3_, 
                               v4_, v5_, F, dyr_, ddyr_, dddyr_);
        return u;
    }

private:
    Controller(const Controller& rhs);
    Controller& operator=(const Controller& rhs);
    
    // sample time
    double dt_;
    // adaptive parameter
    double theta0_;
    double theta_;
    
    // constant
    const double a1_ = 1/m;
    const double a2_ = B/m;
    const double a3_ = 4*beta_e*Kq*Ka/V1;
    const double a4_ = 4*beta_e*A1/V1;
    const double c1_ = 200;
    const double c2_ = 200;
    const double c3_ = 200;
    const double K1_ = -c1_;  // d_alpha1/d_x1
    const double K2_ = c1_;  // d_alpha1/d_yr
    const double K3_ = 1;
    const double v1_ = -(c1_ * c2_ + 1) / a1_;  // d_alpha2/d_x1
    const double v2_ = -(c1_ - a2_ + c2_) / a1_; // d_alpha2/d_x2
    const double v3_ = (c1_ * c2_ + 1) / a1_;  // d_alpha2/d_yr
    const double v4_ = (c2_ + c1_) / a1_;  // d_alpha2/d_dyr
    const double v5_ = 1/a1_;  // d_alpha2/d_ddyr
    const double k_ = 0.05;  // sigma
    const double gamma_ = 0.1;

    // other vars
    double yr_ = 0;
    double dyr_ = 0;
    double ddyr_ = 0;
    double dddyr_ = 0;
    double last_x1_ = 0;

    void updateTheta(double z3, double x3, double theta) {
        double dtheta = gamma_ * (z3*x3 - k_*(theta - theta0_));
        theta_ = theta + dtheta * dt_;
    }

    double calcX2(double x1) {
        double x2 = (x1-last_x1_)/dt_;
        last_x1_ = x1;  // update
        return x2;
    }

    void updateYrDeri(double yr) {
        double last_yr = yr_;
        double last_dyr = dyr_;
        double last_ddyr = ddyr_;
        dyr_ = (yr-last_yr) / dt_;
        ddyr_ = (dyr_ - last_dyr) / dt_;
        dddyr_ = (ddyr_ - last_ddyr) / dt_;
        yr_ = yr;
    }

    double calcZ1(double x1, double yr) const {
        return x1 - yr;
    }
        
    double calcZ2(double x2, double alpha1) const {
        return x2 - alpha1;
    }

    double calcZ3(double x3, double alpha2) const {
        return x3 - alpha2;
    }

    double clacAlpha1(double z1, double dyr) const {
        double alpha1 = -c1_ * z1 + dyr;
        return alpha1;
    }

    double calcAlpha2(double z1, double z2, double x2, double F, double dyr, double ddyr) const {
        double alpha2 = (1/a1_) * (-z1 - c2_*z2 + a2_*x2+F +
                                    K1_*x2 + K2_*dyr + K2_*ddyr);
        return alpha2;
    }

    double calcU(double z2, double z3, double x2, double theta, double x3, double v1, double v2, 
                 double v3, double v4, double v5, double F, double dyr, double ddyr, double dddyr) const {
        double u = (1/a3_) * (-a1_*z2-c3_*z3 + a4_*x2 + theta*x3 + v1*x2 + v2*a1_*x3 -
                              a2_*v2*x2 - v2*F + v3*dyr + v4*ddyr + v5*dddyr);
        return u;
    }

};

#endif