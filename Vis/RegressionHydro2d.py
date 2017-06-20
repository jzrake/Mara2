

class MethodOfLinesTVD(object):

    def plot(self):
        figs = [ ]
        for scheme in ['pcm2', 'plm2']:
            try:
                figs.append(self.with_scheme(scheme))
            except IOError as e:
                print e, "for scheme", scheme
        return figs


class DensityWave2D(MethodOfLinesTVD):

    def with_scheme(self, scheme):
        import os
        import h5py
        import matplotlib.pyplot as plt

        base = 'DensityWave2D-{0}'.format(scheme)
        chkpt0 = h5py.File (os.path.join('data', '{0}.0000.h5'.format(base)), 'r')
        chkpt1 = h5py.File (os.path.join('data', '{0}.0001.h5'.format(base)), 'r')
        d0 = chkpt0['primitive']['density'][...]
        d1 = chkpt1['primitive']['density'][...]

        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        cax = ax1.imshow(d1, interpolation='nearest')
        fig.colorbar(cax)
        fig.suptitle(base)
        return fig


def search(name, terms):
    if not terms: return True

    for term in terms:
        if term.startswith('~') and term[1:] in name:
            return False
        if not term.startswith('~') and term in name:
            return True
    return False


def main():
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("search_terms", nargs='*', help="search terms to match in test figures")
    parser.add_argument("-l", "--list", action='store_true', help="only show the list of tests")
    parser.add_argument("--pdf", action='store_true', help="export figures to PDF")
    args = parser.parse_args()

    plotters = [DensityWave2D()]

    for plotter in plotters:

        name = plotter.__class__.__name__

        if not search(name, args.search_terms):
            continue

        print plotter.__class__

        if args.list:
            continue

        try:
            figs = plotter.plot()

            if args.pdf:
                [fig.savefig(name + '.pdf') for fig in figs]

        except IOError as e:
            print e

    if not args.pdf:
        plt.show()

if __name__ == "__main__":
    main()
