#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Generate plot of the Betamodel for WVTICs paper """

import numpy
import astropy
import matplotlib
matplotlib.rcParams.update({ "figure.figsize": (12, 9), "font.size": 28, })
from matplotlib import pyplot

from parse_gadget_format_two import eat_snapshot

__author__ = "Timo Halbesma"
__email__ = "halbesma@MPA-Garching.MPG.DE"


def p2(a):
    return ( (a) * (a) )


def betamodel(r, rho0=1e-26, beta=2.0/3.0, rc=20):
    return rho0 * pow(1 + p2(r/rc), -3.0/2.0*beta);


def plot_galaxy_cluster():
    pyplot.figure(figsize=(12, 9))
    save_as = "out/IC_GalaxyCluster.png"

    for fname, color, label in zip(
            ["IC_GalaxyCluster_001", "IC_GalaxyCluster_008", "IC_GalaxyCluster_016", "IC_GalaxyCluster_054"],
            ["black", "blue", "brown", "red"], ["Rejection", "8 steps", "16 steps", "54 steps"]):
        header, gas, dm = eat_snapshot(fname)

        boxsize = header["boxSize"]
        boxhalf = boxsize / 2
        gas["r"] = numpy.sqrt( p2(gas["x"] - boxhalf) + p2(gas["y"] - boxhalf) + p2(gas["z"] - boxhalf) )
        pyplot.scatter(gas["r"], gas["rho"], color=color, s=1, facecolor="none", label=label)

    rsample = 1000
    r = numpy.power(10, numpy.linspace(numpy.log10(1), numpy.log10(rsample), 42))
    pyplot.plot(r, betamodel(r), "-", c="orange", label="Model")

    pyplot.xlim(1, boxsize)
    pyplot.ylim(3e-30, 2e-26)
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.xlabel("Radius [kpc]")
    pyplot.ylabel(r"Density [g/cm$^3$]")

    pyplot.tight_layout()
    pyplot.legend(loc="lower left")
    pyplot.savefig(save_as, dpi=300)
    pyplot.close()


if __name__ == "__main__":
    plot_galaxy_cluster()
