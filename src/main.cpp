#include <vector>
#include <cmath>
#include "model.hpp"
#include "observer.hpp"
#include "controller.hpp"


void print_vector(const std::vector<double> &output_data) {
    for(auto ele : output_data) {
        std::cout<<ele<<"\t"; 
    }
    std::cout<<std::endl;
}


int main() {
    const double dt = 0.0001;
    const double sim_time = 10;
    const size_t sim_step = static_cast<size_t>(sim_time/dt);

    Observer observer(dt);
    Model model(dt, 0, 0, 0);  // x1_0=0, x2_0=0, x3_0=0
    Controller controller(8.5, dt); // theta0 = 8.5

    // prepare input data
    std::vector<double> fake_f; // fake load
    std::vector<double> yr_arr;  // fake input
    std::vector<double> dyr_arr;
    std::vector<double> ddyr_arr;
    std::vector<double> dddyr_arr;
    for(size_t i=0; i<sim_step; i++) {
        double t = i*dt;
        double f = std::sin(0.01*t) * 100/2.2;
        fake_f.push_back(f);

        double yr = 0.05 * std::sin(10 * M_PI * t);
        yr_arr.push_back(yr);

        double dyr = 0.05 * 10 * M_PI * std::cos(10 * M_PI * t);
        dyr_arr.push_back(dyr);

        double ddyr = -0.05 * (10 * M_PI) * (10 * M_PI) * std::sin(10 * M_PI * t);
        ddyr_arr.push_back(ddyr);

        double dddyr = -0.05 * (10 * M_PI) * (10 * M_PI) * (10 * M_PI) * std::cos(10 * M_PI * t);
        dddyr_arr.push_back(dddyr);
    }

    // run controller
    std::vector<double> output_arr;
    for(size_t i=0; i<sim_step; i++) {
        // fake sensor
        double x1 = model.getX1();
        double x2 = model.getX2();
        double x3 = model.getX3();
        // update observer
        observer.update(x1, x3);
        double est_f = observer.getF();

        // run controller
        double u = controller.run(yr_arr[i], x1, x2, x3, est_f, dyr_arr[i], ddyr_arr[i], dddyr_arr[i]);
        output_arr.push_back(u);  // store output

        // update
        model.update(fake_f[i], u);
    }

    std::cout<<"output data:" << std::endl;
    std::cout<<output_arr.back() << std::endl;
    // print_vector(output_arr);
}