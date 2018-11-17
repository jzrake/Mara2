import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter



def rtheta_mesh(fig, filename):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    h5f = h5py.File(filename, 'r')

    r = h5f['mesh']['points']['x'][:]
    q = h5f['mesh']['points']['y'][:]
    p  = h5f['primitive']['pressure'][:].T
    d  = h5f['primitive']['density'][:].T
    v1 = h5f['primitive']['velocity1'][:].T
    v2 = h5f['primitive']['velocity2'][:].T

    r /= r[0]
    u1 = v1 / (1 - v1 * v1 - v2 * v2)**0.5

    R, Q = np.meshgrid(r, q)
    X = np.log10(R) * np.cos(Q)
    Y = np.log10(R) * np.sin(Q)

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    div1 = make_axes_locatable(ax1)
    div2 = make_axes_locatable(ax2)
    cax1 = div1.append_axes('right', size='5%', pad=0.05)
    cax2 = div2.append_axes('right', size='5%', pad=0.05)

    im1 = ax1.pcolormesh(Y, X, np.log10(d), edgecolor='none')
    im2 = ax2.pcolormesh(Y, X, u1, edgecolor='none')

    fig.suptitle(r"$t = {:0.2f} ms$".format(1e3 * h5f['status']['simulationTime'].value))
    fig.colorbar(im1, cax=cax1, orientation='vertical')
    fig.colorbar(im2, cax=cax2, orientation='vertical')
    fig.subplots_adjust(left=0.08, right=0.96, top=1.0, bottom=0.00, wspace=0.25)

    ax1.set_aspect('equal')
    ax1.set_title(r"$\log_{10} \rho$")
    ax1.set_xlabel(r'$\log_{10} (x/50 \rm{km})$')
    ax1.set_ylabel(r'$\log_{10} (y/50 \rm{km})$')

    ax2.set_aspect('equal')
    ax2.set_title(r"$u_r$")
    ax2.set_xlabel(r'$\log_{10} (x/50 \rm{km})$')
    ax2.set_ylabel(r'$\log_{10} (y/50 \rm{km})$')



def rtheta_mesh_plot(args):
    fig = plt.figure(figsize=[10,6])
    rtheta_mesh(fig, args.filenames[0])
    plt.show()



def rtheta_mesh_movie(args):

    writer = FFMpegWriter(fps=15)
    fig = plt.figure(figsize=[10,6])
    dpi = 200

    with writer.saving(fig, "out.mp4", dpi):
        for filename in args.filenames:
            print(filename)
            rtheta_mesh(fig, filename)
            writer.grab_frame()
            fig.clf()



def radial_profile(fig, filename):
    h5f = h5py.File(filename, 'r')

    r = h5f['mesh']['points']['x'][:]
    r = 0.5 * (r[1:] + r[:-1])
    vr = h5f['primitive']['velocity1'][:]
    dg = h5f['primitive']['density'][:]
    pg = h5f['primitive']['pressure'][:]

    fig.suptitle(r"$t = {:0.1f} \ r_0/c$".format(h5f['status']['simulationTime'].value))
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)
    ax1.plot(r, vr / (1 - vr**2)**0.5, label=r'$u_r$')
    ax2.plot(r, dg, label=r'$\rho$')
    ax3.plot(r, pg, label=r'$p$')

    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax3.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    # ax1.set_ylim(5e-2, 1e4)
    # ax2.set_ylim(1e-24, 10)
    # ax3.set_ylim(1e-24, 10)
    ax1.set_ylabel(r'four-velocity $u_r$')
    ax2.set_ylabel(r'density $\rho$')
    # ax2.set_xlabel(r'$r / r_0$')
    ax3.set_ylabel(r'pressure $p$')
    ax3.set_xlabel(r'$r / r_0$')



def radial_profile_plot(args):

    fig = plt.figure(figsize=[6,8])

    radial_profile(fig, args.filenames[0])
    plt.show()



def radial_profile_movie(args):
    import os

    fig = plt.figure(figsize=[6,8])

    for filename in args.filenames:
        pngname = filename.replace('.h5', '.png')
        print(pngname)
        radial_profile(fig, filename)
        fig.savefig(pngname)
        fig.clf()

    os.system("ffmpeg -i jet-400/chkpt.%04d.png -f mp4 -vcodec h264 -pix_fmt yuv420p out.mp4")



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['radial-plot', 'radial-movie', 'rtheta-plot', 'rtheta-movie'])
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('--equator', action='store_true')
    parser.add_argument('--fixedr', action='store_true')
    args = parser.parse_args()

    commands = {
    'radial-plot':radial_profile_plot, 'radial-movie':radial_profile_movie,
    'rtheta-plot':rtheta_mesh_plot, 'rtheta-movie':rtheta_mesh_movie
    }

    commands[args.command](args)
