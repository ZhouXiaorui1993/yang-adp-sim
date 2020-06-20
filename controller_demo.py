#!/usr/bin/env/python3
# -*- coding: utf-8 -*-
from  config import *
from Model_demo import Model, Observer


class Controller:
    def __init__(self, theta0, dt):
        # param init
        # sample time
        self.dt = dt

        # adaptive params
        self.theta0 = theta0
        self.theta = theta0

        self.a1 = 1/m
        self.a2 = B/m
        self.a3 = 4*beta_e*Kq*Ka/V1
        self.a4 = 4*beta_e*A1/V1

        self.c1 = 200
        self.c2 = 200
        self.c3 = 200
        self.K1 = -self.c1  # d_alpha1/d_x1
        self.K2 = self.c1  # d_alpha1/d_yr
        self.K3 = 1
        self.v1 = -(self.c1 * self.c2 + 1) / self.a1  # d_alpha2/d_x1
        self.v2 = -(self.c1 - self.a2 + self.c2) / self.a1  # d_alpha2/d_x2
        self.v3 = (self.c1 * self.c2 + 1) / self.a1  # d_alpha2/d_yr
        self.v4 = (self.c2 + self.c1) / self.a1  # d_alpha2/d_dyr
        self.v5 = 1/self.a1  # d_alpha2/d_ddyr

        self.k = 0.05  # sigma
        self.gamma = 0.1

        # other vars
        self.yr = 0
        self.dyr = 0
        self.ddyr = 0
        self.dddyr = 0
        self.last_x1 = 0
        # self.epsilon = 1

    def _update_theta(self, z3, x3, theta):
        dtheta = self.gamma * (z3*x3 - self.k*(theta - self.theta0))
        self.theta = self.theta + dtheta*self.dt

    def _calc_x2(self, x1):
        x2 = (x1-self.last_x1)/self.dt
        self.last_x1 = x1  # update
        return x2

    def _update_yr_deri(self, yr):
        last_yr = self.yr
        last_dyr = self.dyr
        last_ddyr = self.ddyr
        self.dyr = (yr-last_yr)/self.dt
        self.ddyr = (self.dyr-last_dyr)/self.dt
        self.dddyr = (self.ddyr-last_ddyr)/self.dt
        self.yr = yr

    @staticmethod
    def _calc_z1(x1, yr):
        return x1 - yr

    @staticmethod
    def _calc_z2(x2, alpha1):
        return x2 - alpha1

    @staticmethod
    def _calc_z3(x3, alpha2):
        return x3 - alpha2

    def _calc_alpha1(self, z1, dyr):
        alpha1 = -self.c1 * z1 + dyr
        return alpha1

    def _calc_alpha2(self, z1, z2, x2, F, dyr, ddyr):
        alpha_2 = (1/self.a1) * (-z1 - self.c2*z2 + self.a2*x2+F +
                                 self.K1*x2 + self.K2*dyr + self.K2*ddyr)

        return alpha_2

    def _calc_u(self, z2, z3, x2, theta, x3, v1, v2, v3, v4,
                v5, F, dyr, ddyr, dddyr):
        u = (1/self.a3) * (-self.a1*z2-self.c3*z3 + self.a4*x2 +
                           theta*x3 + v1*x2 + v2*self.a1*x3 -
                           self.a2*v2*x2 - v2*F + v3*dyr + v4*ddyr + v5*dddyr)

        return u

    def run(self, yr, x1, x2, x3, F, dyr, ddyr, dddyr):
        # update
        # self._update_yr_deri(yr)
        self.dyr = dyr
        self.ddyr = ddyr
        self.dddyr = dddyr

        # calculation
        z1 = self._calc_z1(x1, yr)
        alpha1 = self._calc_alpha1(z1, self.dyr)
        # x2 = self._calc_x2(x1)
        z2 = self._calc_z2(x2, alpha1)
        alpha2 = self._calc_alpha2(z1, z2, x2, F, self.dyr, self.ddyr)
        z3 = self._calc_z3(x3, alpha2)

        self._update_theta(z3, x3, self.theta)

        u = self._calc_u(z2, z3, x2, self.theta, x3, self.v1, self.v2, self.v3, self.v4, self.v5,
                         F, self.dyr, self.ddyr, self.dddyr)
        return u


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    dt = 0.0001
    sim_time = 10
    sim_step = int(sim_time/dt)

    observer = Observer(dt)
    model = Model(dt, x1_0=0, x2_0=0, x3_0=0)  # TODO: set correct init value
    theta0 = 8.5
    control_obj = Controller(theta0, dt)

    # fake load
    t = np.arange(sim_step) * dt
    f = np.sin(0.01*t) * 100/2.2

    # fake input
    yr = 0.05 * np.sin(10 * np.pi * t)
    # debug
    dyr_lst = 0.05 * 10 * np.pi * np.cos(10 * np.pi * t)
    ddyr_lst = -0.05 * (10 * np.pi) * (10 * np.pi) * np.sin(10 * np.pi * t)
    dddyr_lst = -0.05 * (10 * np.pi) * (10 * np.pi) * (10 * np.pi) * np.cos(10 * np.pi * t)
    # yr = [0.05]*sim_step

    output_lst = []

    for i in range(sim_step):
        # fake sensor
        x1, x2, x3 = model.get_state()
        # update observer
        observer.update_f(x1, x3)
        est_f = observer.get_f()
        # run controller
        u = control_obj.run(yr[i], x1, x2, x3, est_f, dyr_lst[i], ddyr_lst[i], dddyr_lst[i])
        # update model
        model.update(f[i], u)
        output_lst.append(model.get_y())

    x = list(range(sim_step))
    assert len(x) == sim_step
    assert len(yr) == sim_step
    assert len(output_lst) == sim_step

    plt.figure()
    plt.plot(x, yr, "-", label="yr")
    plt.plot(x, output_lst, "-.", label="output")
    # plt.plot(x, dyr_lst, label="dyr")
    # plt.plot(x, ddyr_lst, label="ddyr")
    # plt.plot(x, dddyr_lst, label="dddyr")
    plt.legend()
    plt.show()