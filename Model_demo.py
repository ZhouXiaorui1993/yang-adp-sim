#!/usr/bin/env/python3
# -*- coding: utf-8 -*-
from config import *


class Model:
    def __init__(self, dt, x1_0, x2_0, x3_0):
        self.x1 = x1_0
        self.x2 = x2_0
        self.x3 = x3_0
        self.dt = dt

        self.a1 = 1/m
        self.a2 = B/m
        self.a3 = 4*beta_e*Kq*Ka/V1
        self.a4 = 4*beta_e*A1/V1

        self.theta = 4*beta_e*Ct/V1  # a5

    def update(self, f, u):
        """
        update model states based on control cmd
        :param dt: sample time
        :param f:
        :param u:
        :return:
        """
        x1 = self.x1
        x2 = self.x2
        x3 = self.x3
        self.x1 = x1 + x2 * self.dt
        self.x2 = x2 + (self.a1*x3 - self.a2*x2 - f) * self.dt
        self.x3 = x3 + (self.a3*u - self.a4*x2 - self.theta*x3) * self.dt

    def get_state(self):
        return self.x1, self.x2, self.x3

    def get_y(self):
        return self.x1


class Observer:
    def __init__(self, dt):
        self.dt = dt
        # constant
        self.k0 = 1000
        self.lmbd = 0  # lambda
        self.m = 2.2

        # var
        self.f = 0
        self.last_x1 = 0  # TODO:

    def update_f(self, x1, x3):
        p = self.k0 * self.m * (x1 - self.last_x1)/self.dt
        d_lambda = -self.k0 * (self.lmbd + p) + self.k0 * x3
        df = self.lmbd + p
        # update f
        self.f = self.f + df * self.dt
        # update lambda for next use
        self.lmbd = self.lmbd + d_lambda * self.dt
        # store x1 for next use
        self.last_x1 = x1

    def get_f(self):
        return self.f
