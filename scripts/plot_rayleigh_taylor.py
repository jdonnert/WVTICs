#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Generate plots to inspect Rayleigh Taylor Instability ICs """

import h5py
import numpy
import pandas
from matplotlib import pyplot

from parse_gadget_format_two import eat_snapshot

__author__ = "Timo Halbesma"
__email__ = "halbesma@MPA-Garching.MPG.DE"


def p2(a):
    return ( (a) * (a) )


def density(x):
    rho1 = 1.0
    rho2 = 2.0
    delta = 0.025

    return rho1 + ( rho2 - rho1 ) / ( 1 + numpy.exp ( - ( x - 0.5 ) / delta ) )


def velocity_x(x, y):
    vx = numpy.zeros(x.shape)

    mask = numpy.where((x > 0.3) & (x < 0.7))
    vx[mask] = 0.025 * ( 1 + numpy.cos ( 8 * numpy.pi * ( y[mask] + 0.25 ) ) ) \
        * ( 1 + numpy.cos( 5 * numpy.pi * ( x[mask] - 0.5 ) ) )

    return vx


def magnetic_field_y(x):
    return (0.07 * numpy.sqrt ( 4 * numpy.pi )) * numpy.ones(x.shape)


def internal_energy(x):
    gamma = 1.4
    rho1 = 1.0
    rho2 = 2.0
    delta = 0.025
    grav_acc = -0.5

    pressure = rho2 / gamma + grav_acc * density(x) * ( x - 0.5 )

    return pressure / ( gamma - 1 ) / density (x)


def plot_rayleigh_taylor_ics():
    save_as = "out/IC_RayleighTaylor.png"

    fname = "../WVTICs_runs/1E5/RayleighTaylor/IC_RayleighTaylorInstability"
    header, gas, dm = eat_snapshot(fname)
    # Caution, boxSize in header is 1D
    # boxsize = header["boxSize"]
    boxsize = numpy.array([1.0, 0.5, 0.1])
    boxhalf = boxsize / 2

    # To check, we take the ICs from Phil Hopkins, taken from
    # http://www.tapir.caltech.edu/~phopkins/sims/ accessed at 20171006
    fname = "rt_ics.hdf5"  # does not have magnetic field
    hopkins_ics = h5py.File("rt_ics.hdf5",'r+')["PartType0"]
    hopkins_pos = hopkins_ics["Coordinates"]
    hopkins_rho = hopkins_ics["Density"]
    hopkins_vel = hopkins_ics["Velocities"]
    hopkins_uint = hopkins_ics["InternalEnergy"]

    # To check sampled ICs, we compare to analytical requirement
    npart = len(hopkins_pos)
    pos = numpy.array([boxsize[0]*numpy.random.random(npart),
        boxsize[1]*numpy.random.random(npart),
        boxsize[2]*numpy.random.random(npart)]).T
    rho = density(pos[::,0])
    By = magnetic_field_y(pos[::,0])
    vx = velocity_x(pos[::,0], pos[::,1])
    uint = internal_energy(pos[::,0])

    fig, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2, figsize=(16, 12))

    # Plot the density
    ax1.plot(gas["x"], gas["rho"], "kD", ms=1, c="g", label="WVTICs", alpha=0.1)
    ax1.plot(hopkins_pos[::,1], hopkins_rho, "kD", ms=1, c="k", label="Hopkins")
    ax1.plot(pos[::,0], rho, "kD", ms=0.5, c="r", alpha=0.1, label="Required")

    # Hack legend to show solid line rather than an unreadable dot
    handles, labels = ax1.get_legend_handles_labels()
    wvticsArtist = pyplot.Line2D((0,1), (0,0), c="g", ls="-")
    hopkinsArtist = pyplot.Line2D((0,1), (0,0), c="k", ls="-")
    requiredArtist = pyplot.Line2D((0,1), (0,0), c="r", ls="-")
    ax1.legend([wvticsArtist, hopkinsArtist, requiredArtist], labels, loc="upper left")

    # Plot the magnetic field in y-direction
    ax2.plot(gas["x"], gas["By"], "kD", ms=1, c="g", label="WVTICs", alpha=0.1)
    ax2.plot(pos[::,0], By, "kD", ms=0.5, c="r", alpha=0.1, label="Required")
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend([wvticsArtist, requiredArtist], labels, loc="upper right")
    ax2.set_ylim(0, 1)

    # Plot the velocity perturbation
    ax3.plot(gas["x"], gas["vx"], "kD", ms=1, c="g", label="WVTICs", alpha=0.1)
    ax3.plot(hopkins_pos[::,1], hopkins_vel[::,1], "kD", ms=1, c="k", label="Hopkins")
    ax3.plot(pos[::,0], vx, "kD", ms=0.5, c="r", alpha=0.1, label="Required")

    # Hack legend to show solid line rather than an unreadable dot
    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend([wvticsArtist, hopkinsArtist, requiredArtist], labels, loc="upper left")

    # Plot the internal energy
    ax4.plot(gas["x"], gas["u"], "kD", ms=1, c="g", label="WVTICs", alpha=0.1)
    ax4.plot(hopkins_pos[::,1], hopkins_uint, "kD", ms=1, c="k", label="Hopkins")
    ax4.plot(pos[::,0], uint, "kD", ms=0.5, c="r", alpha=0.1, label="Required")

    # Hack legend to show solid line rather than an unreadable dot
    handles, labels = ax4.get_legend_handles_labels()
    ax4.legend([wvticsArtist, hopkinsArtist, requiredArtist], labels, loc="upper right")

    # Set labels and limits
    ax3.set_xlabel("x [WVTICs] == y [Hopkins]")
    ax4.set_xlabel("x [WVTICs] == y [Hopkins]")
    ax4.set_xlim(0, 1)

    ax1.set_ylabel(r"$\rho$")
    ax2.set_ylabel(r"$B_x$")
    ax3.set_ylabel(r"$v_x$")
    ax4.set_ylabel(r"$U$")

    pyplot.tight_layout()
    pyplot.savefig(save_as)
    pyplot.close()


def plot_rayleigh_taylor_diagnostics():
    save_as = "out/IC_RayleighTaylor_diagnostics.png"
    logfile = "../WVTICs_runs/1E5/RayleighTaylor/diagnostics.log"

    data = pandas.read_csv(logfile, quotechar='"',
        skipinitialspace=True, delimiter='\t').as_matrix()

    pyplot.figure()
    ax = pyplot.gca()

    # ax.plot(data[::,0], data[::,1], label="min")
    ax.plot(data[::,0], data[::,2], label="max",
        c="black", alpha=1, lw=1)
    ax.plot(data[::,0], data[::,3], label="mean",
        c="blue", alpha=1, lw=1)
    ax.set_xlabel("#Iterations", fontsize=6)
    ax.set_xlim(0, 500)
    ax.set_ylabel("Error", fontsize=6)
    ax.set_yscale("log")
    ax.set_ylim(0.005, 1.05)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)
        tick.label.set_rotation("vertical")

    pyplot.legend(loc="center right")
    pyplot.tight_layout()
    pyplot.savefig(save_as)
    pyplot.close()


if __name__ == "__main__":
    plot_rayleigh_taylor_ics()
    plot_rayleigh_taylor_diagnostics()
