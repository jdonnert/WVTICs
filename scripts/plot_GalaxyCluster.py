#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Generate plot of the Betamodel for WVTICs paper """

import numpy
import astropy
from matplotlib import pyplot

from parse_gadget_format_two import eat_snapshot

__author__ = "Timo Halbesma"
__email__ = "halbesma@MPA-Garching.MPG.DE"


def p2(a):
    return ( (a) * (a) )


def betamodel(r, rho0=1e-26, beta=2.0/3.0, rc=20):
    return rho0 * pow(1 + p2(r/rc), -3.0/2.0*beta);


def plot_galaxy_cluster(fname, save_as):
    header, gas, dm = eat_snapshot(fname)

    boxsize = header["boxSize"]
    boxhalf = boxsize / 2
    gas["r"] = numpy.sqrt( p2(gas["x"] - boxhalf) + p2(gas["y"] - boxhalf) + p2(gas["z"] - boxhalf) )

    pyplot.figure(figsize=(12, 9))

    pyplot.plot(gas["r"], gas["rho"], "ko", ms=1, markerfacecolor="none")
    rsample = 100
    r = numpy.power(10, numpy.linspace(numpy.log10(1), numpy.log10(rsample), 42))
    pyplot.plot(r, betamodel(r), "r-")

    pyplot.xlim(1, boxsize)
    pyplot.xscale("log")
    pyplot.yscale("log")

    pyplot.tight_layout()
    pyplot.savefig(save_as, dpi=300)
    pyplot.close()


if __name__ == "__main__":
    for fname in ["IC_GalaxyCluster"]:
        plot_galaxy_cluster(fname, "out/IC_GalaxyCluster.png")
