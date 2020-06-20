#ifndef CONTROL_MODEL_H_
#define CONTROL_MODEL_H_

#include "config.h"

class Model {
public:
    Model(double dt, double x1_0, double x2_0, double x3_0) : 
        dt_(dt), x1_(x1_0), x2_(x2_0), x3_(x3_0) {}
    ~Model() {}

    void update(double f, double u) {
        double x1 = x1_;
        double x2 = x2_;
        double x3 = x3_;
        x1_ = x1 + x2 * dt_;
        x2_ = x2 + (a1_ * x3 - a2_ * x2 - f) * dt_;
        x3 = x3 + (a3_ * u - a4_ * x2 - theta_ * x3) * dt_;
    }

    double getX1() const {
        return x1_;
    }

    double getX2() const {
        return x2_;
    }

    double getX3() const {
        return x3_;
    }

    double getY() const {
        // output: y=x1
        return x1_;
    }

private:
    Model(const Model& rhs);
    Model& operator=(const Model& rhs);
    
    // model states
    double  x1_;
    double  x2_;
    double  x3_;
    double  dt_;

    // constant
    const double a1_ = 1 / m;
    const double a2_ = B / m;
    const double a3_ = 4 * beta_e * Kq * Ka / V1;
    const double a4_ = 4 * beta_e * A1 / V1;
    const double theta_ = 4 * beta_e * Ct / V1;  // a5
};

#endif