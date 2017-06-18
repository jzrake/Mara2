

def get_data(checkpoint, field):
    import h5py
    h5f = h5py.File(checkpoint, 'r')

    if field == 'B':
        B1 = h5f['primitive']['magnetic1'][...]
        B2 = h5f['primitive']['magnetic2'][...]
        B3 = h5f['primitive']['magnetic3'][...]
        return (B1**2 + B2**2 + B3**2)**0.5

    elif field == 'd': return h5f['primitive']['density'][...]
    elif field == 'p': return h5f['primitive']['pressure'][...]
    elif field == 'u': return h5f['primitive']['velocity1'][...]
    elif field == 'v': return h5f['primitive']['velocity2'][...]
    elif field == 'w': return h5f['primitive']['velocity3'][...]


def plot_linear(fig, checkpoint, args):
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    y = get_data(checkpoint, args.field)
    ax1.plot(y, '-o', mfc='none')


def plot_image(fig, checkpoint, args):
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    img = get_data(checkpoint, args.field)
    iax = ax1.imshow(img, interpolation='none')
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    fig.colorbar(iax)


def plot_hydro1d(fig, checkpoint, args):
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    d = get_data(checkpoint, 'd')
    p = get_data(checkpoint, 'p')
    u = get_data(checkpoint, 'u')
    ax1.plot(d, '-o', mfc='none', label=r'$\rho$')
    ax1.plot(p, '-o', mfc='none', label=r'$p$')
    ax1.plot(u, '-o', mfc='none', label=r'$u$')
    ax1.legend(loc='best')
    fig.suptitle(checkpoint)


def plot(args):
    import matplotlib.pyplot as plt

    plot_function = {
    'linear': plot_linear,
    'image': plot_image,
    'hydro1d': plot_hydro1d }[args.command]

    for checkpoint in args.checkpoints[::-1]:
        fig = plt.figure()
        plot_function(fig, checkpoint, args)

    plt.show()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('checkpoints', nargs='+')
    parser.add_argument('--field', '-f', default='p', choices=['B', 'v', 'd'])
    args = parser.parse_args()
    plot(args)


if __name__ == "__main__":
    main()
