#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.axes_grid1 import make_axes_locatable



def make_binary_plot(fig, filename):
    h5f = h5py.File(filename, 'r')

    R  = h5f['user']['DomainRadius'][...]
    x  = h5f['mesh']['points']['x'][...]
    y  = h5f['mesh']['points']['y'][...]
    p  = h5f['primitive']['pressure'][...]
    d  = h5f['primitive']['sigma'][...]
    vx = h5f['primitive']['velocity1'][...]
    vy = h5f['primitive']['velocity2'][...]

    xc = 0.5 * (x[1:] + x[:-1])[:,None]
    yc = 0.5 * (y[1:] + y[:-1])[None,:]

    r2 = xc**2 + yc**2
    xh = xc / r2**0.5
    yh = yc / r2**0.5
    vr = vx * xh + vy * yh
    vq = vy * xh - vx * yh
    cs = (p / d)**0.5
    # Mq = vq / cs

    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    div1 = make_axes_locatable(ax1)
    div2 = make_axes_locatable(ax2)
    div3 = make_axes_locatable(ax3)
    div4 = make_axes_locatable(ax4)

    cax1 = div1.append_axes('bottom', size='5%', pad=0.05)
    cax2 = div2.append_axes('bottom', size='5%', pad=0.05)
    cax3 = div3.append_axes('bottom', size='5%', pad=0.05)
    cax4 = div4.append_axes('bottom', size='5%', pad=0.05)

    im1 = ax1.imshow(d,  extent=[-R, R, -R, R])#, vmin=0, vmax=1)
    im2 = ax2.imshow(p,  extent=[-R, R, -R, R])#, vmin=0, vmax=0.001)
    im3 = ax3.imshow(vq, extent=[-R, R, -R, R])#, vmin=-0.5, vmax=1.5)
    im4 = ax4.imshow(vr, extent=[-R, R, -R, R])#, vmin=-0.5, vmax=1.5)

    fig.colorbar(im1, cax=cax1, orientation='horizontal')
    fig.colorbar(im2, cax=cax2, orientation='horizontal')
    fig.colorbar(im3, cax=cax3, orientation='horizontal')
    fig.colorbar(im4, cax=cax4, orientation='horizontal')

    ax1.set_title(r"$\rho$")
    ax2.set_title(r"$P$")
    # ax3.set_title(r"$\mathcal{M} = v_\theta / c_s$")
    ax3.set_title(r"$v_\theta$")
    ax4.set_title(r"$v_r$")

    t = np.linspace(0, 2 * np.pi, 300)
    x = 7. * np.cos(t)
    y = 7. * np.sin(t)
    ax4.plot(x, y, '--', lw=1, c='k')

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(r"$t = {:.1f} \ t_{{\rm orb}}$".format(h5f['status']['simulationTime'].value / (2 * np.pi)))
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.15, wspace=0.15)



def binary_plot(args):
    for filename in args.filenames:
        fig = plt.figure(figsize=[8,8])
        make_binary_plot(fig, filename)
    plt.show()



def binary_movie(args):

    writer = FFMpegWriter(fps=15)
    fig = plt.figure(figsize=[8,8])
    dpi = 200

    with writer.saving(fig, args.output, dpi):
        for filename in args.filenames:
            print(filename)
            make_binary_plot(fig, filename)
            writer.grab_frame()
            fig.clf()



def make_shear_plot1d(ax, filename):
    h5f = h5py.File(filename, 'r')

    nu = 0.1
    t0 = 1.0 / nu
    t  = h5f['status']['simulationTime'].value + t0
    x  = h5f['mesh']['points']['x'][...]
    y  = h5f['mesh']['points']['y'][...]
    p  = h5f['primitive']['pressure'][...]
    d  = h5f['primitive']['density'][...]
    vx = h5f['primitive']['velocity1'][...]
    vy = h5f['primitive']['velocity2'][...]

    xc = 0.5 * (x[1:] + x[:-1])
    f = (t0 / t)**0.5 * np.exp(-xc**2 / (4 * nu * t))

    ax.plot(xc, vy, 'o', mfc='none', label=filename)
    ax.plot(xc, f, '-', mfc='none')
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$u_y$")



def shear_plot1d(args):
    fig = plt.figure(figsize=[8,6])
    ax1 = fig.add_subplot(1, 1, 1)

    for filename in args.filenames:
        make_shear_plot1d(ax1, filename)

    fig.legend()
    plt.show()



def make_shear_plot2d(ax1, ax2, filename):
    h5f = h5py.File(filename, 'r')

    nu = 0.1
    t0 = 1.0 / nu
    t  = h5f['status']['simulationTime'].value + t0
    x  = h5f['mesh']['points']['x'][...]
    y  = h5f['mesh']['points']['y'][...]
    p  = h5f['primitive']['pressure'][...]
    d  = h5f['primitive']['density'][...]
    vx = h5f['primitive']['velocity1'][...]
    vy = h5f['primitive']['velocity2'][...]

    xc = 0.5 * (x[1:] + x[:-1])[:,None]
    yc = 0.5 * (y[1:] + y[:-1])[None,:]
    q = 0.0
    nx = np.cos(q)
    ny = np.sin(q)
    r = xc * nx + yc * ny
    f = (t0 / t)**0.5 * np.exp(-r**2 / (4 * nu * t))

    ax1.imshow(vx)
    ax2.imshow(vy)
    ax1.set_title(r"$u_x$")
    ax2.set_title(r"$u_y$")

    div1 = make_axes_locatable(ax1)
    div2 = make_axes_locatable(ax2)

    cax1 = div1.append_axes('bottom', size='5%', pad=0.05)
    cax2 = div2.append_axes('bottom', size='5%', pad=0.05)

    im1 = ax1.imshow(vx.T)
    im2 = ax2.imshow(vy.T)

    for ax in [ax1, ax2]:
        ax.set_xticks([])
        ax.set_yticks([])

    return im1, im2, cax1, cax2



def shear_plot2d(args):
    for filename in args.filenames:
        fig = plt.figure(figsize=[8,6])
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
    
        im1, im2, cax1, cax2 = make_shear_plot2d(ax1, ax2, filename)
    
        fig.suptitle(filename)
        fig.colorbar(im1, cax=cax1, orientation='horizontal')
        fig.colorbar(im2, cax=cax2, orientation='horizontal')
    plt.show()



def make_radial_plot(ax, filename):
    h5f = h5py.File(filename, "r")

    x  = h5f['mesh']['points']['x'][...]
    y  = h5f['mesh']['points']['y'][...]
    p  = h5f['primitive']['pressure'][...]
    d  = h5f['primitive']['density'][...]
    vx = h5f['primitive']['velocity1'][...]
    vy = h5f['primitive']['velocity2'][...]

    xc = 0.5 * (x[1:] + x[:-1])[:,None]
    yc = 0.5 * (y[1:] + y[:-1])[None,:]

    r2 = xc**2 + yc**2
    xh = xc / r2**0.5
    yh = yc / r2**0.5
    vr = vx * xh + vy * yh
    vq = vy * xh - vx * yh

    N = x.shape[0] - 1

    ax.plot(yc[0, N//2:], d[N // 2, N // 2:])



def radial_plot(args):
    fig = plt.figure(figsize=[8,6])
    ax = fig.add_subplot(1, 1, 1)

    for filename in args.filenames:
        make_radial_plot(ax, filename)

    plt.show()



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['binary-plot', 'binary-movie', 'shear-plot1d', 'shear-plot2d', 'radial-plot'])
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('--output', '-o', default='out.mp4')
    args = parser.parse_args()

    commands = {
    'binary-plot':binary_plot,
    'binary-movie':binary_movie,
    'shear-plot1d':shear_plot1d,
    'shear-plot2d':shear_plot2d,
    'radial-plot':radial_plot
    }

    commands[args.command](args)
