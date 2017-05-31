
def get_image_data(checkpoint, field):
    import h5py
    h5f = h5py.File(checkpoint, 'r')

    if field == 'B':
        B1 = h5f['primitive']['magnetic1'][...]
        B2 = h5f['primitive']['magnetic2'][...]
        B3 = h5f['primitive']['magnetic3'][...]
        return (B1**2 + B2**2 + B3**2)**0.5

    elif field == 'd':
        return h5f['primitive']['density']

def plot(args):
    import matplotlib.pyplot as plt

    for checkpoint in args.checkpoints:
        fig = plt.figure()
        ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        img = get_image_data(checkpoint, args.field)
        iax = ax1.imshow(img, interpolation='none')
        ax1.xaxis.set_ticklabels([])
        ax1.yaxis.set_ticklabels([])
        fig.colorbar(iax)

    plt.show()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('checkpoints', nargs='+')
    parser.add_argument('--field', '-f', default='B', choices=['B', 'v', 'd'])
    args = parser.parse_args()
    plot(args)


if __name__ == "__main__":
    main()
