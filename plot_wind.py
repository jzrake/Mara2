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

    ax = fig.add_subplot(1, 1, 1)
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='5%', pad=0.05)

    #im = ax.pcolormesh(Y, X, np.log10(d), edgecolor='none')#, vmin=-1.0, vmax=1.4)

    im = ax.pcolormesh(Y, X, u1, edgecolor='none')#, vmin=-1.0, vmax=1.4)

    alpha = h5f['user']['alpha'].value
    theta0 = h5f['user']['theta0'].value

    q = np.linspace(0, np.pi / 2, 128)
    R = np.exp(-(q / theta0)**alpha)
    ax.plot(R * np.sin(q), R * np.cos(q), c='k', lw=2, ls='--')

    fig.suptitle(r"$\log_{{10}} u_r$, $t = {:0.1f} \ r_0/c$".format(h5f['status']['simulationTime'].value))
    fig.colorbar(im, cax=cax, orientation='vertical')
    fig.subplots_adjust(top=0.98, bottom=0.02)

    ax.set_aspect('equal')
    ax.set_xlabel(r'$\log_{10} x / r_0$')
    ax.set_ylabel(r'$\log_{10} y / r_0$')



def rtheta_mesh_plot(args):
    fig = plt.figure(figsize=[8,8])
    rtheta_mesh(fig, args.filenames[0])
    plt.show()



def rtheta_mesh_movie(args):

    writer = FFMpegWriter(fps=15)
    fig = plt.figure(figsize=[8,8])
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
