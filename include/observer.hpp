#ifndef CONTROL_OBSERVER_H_
#define CONTROL_OBSERVER_H_

#include "config.h"


class Observer{
public:
    Observer(double dt): dt_(dt) {}
    ~Observer() {}

    double getF() const {
        return f_;
    }

    void update(double x1, double x3) {
        const double p = -k0_ * m * (x1 - last_x1_) / dt_;
        const double d_lambda = -k0_ * (lmbd_ + p) + k0_ * x3;

        const double df = lmbd_ + p;
        // update f
        f_ = f_ + df * dt_;
        // update lambda for next use
        lmbd_ = lmbd_ + d_lambda * dt_;
        // store x1 for next use
        last_x1_ = x1;
    }

private:
    Observer(const Observer& rhs);
    Observer& operator=(const Observer& rhs);

    double dt_;
    // constant
    const double k0_ = 1000;

    // var
    double f_ = 0;  // output
    double last_x1_ = 0;
    double lmbd_ = 0; // lambda
};

#endif