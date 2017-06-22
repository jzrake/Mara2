

class MethodOfLinesTVD(object):

    def plot(self):
        figs = [ ]
        for scheme in ['pcm1', 'pcm2', 'pcm3', 'plm1', 'plm2', 'plm3']:
            try:
                figs.append(self.with_scheme(scheme))
            except IOError as e:
                print e
        return figs


class Shocktube1(MethodOfLinesTVD):

    def with_scheme(self, scheme):
        import os
        import h5py
        import matplotlib.pyplot as plt

        base = 'Shocktube1-{0}'.format(scheme)
        chkpt = h5py.File (os.path.join('data', '{0}.0001.h5'.format(base)), 'r')
        d = chkpt['primitive']['density'][...]
        p = chkpt['primitive']['pressure'][...]
        u = chkpt['primitive']['velocity1'][...]

        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax1.plot(d, '-o', mfc='none', label=r'$\rho$')
        ax1.plot(p, '-o', mfc='none', label=r'$p$')
        ax1.plot(u, '-o', mfc='none', label=r'$u$')
        ax1.legend(loc='best')
        fig.suptitle(base)
        return fig


class Shocktube3(MethodOfLinesTVD):

    def with_scheme(self, scheme):
        import os
        import h5py
        import matplotlib.pyplot as plt

        base = 'Shocktube3-{0}'.format(scheme)
        chkpt = h5py.File (os.path.join('data', '{0}.0001.h5'.format(base)), 'r')
        d = chkpt['primitive']['density'][...]
        p = chkpt['primitive']['pressure'][...]
        u = chkpt['primitive']['velocity1'][...]

        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax1.plot(d, '-o', mfc='none', label=r'$\rho$')
        ax1.plot(p, '-o', mfc='none', label=r'$p$')
        ax1.set_yscale('log')
        ax1.legend(loc='best')
        fig.suptitle(base)
        return fig


class DensityWave(MethodOfLinesTVD):

    def with_scheme(self, scheme):
        import os
        import h5py
        import matplotlib.pyplot as plt

        base = 'DensityWave-{0}'.format(scheme)
        chkpt0 = h5py.File (os.path.join('data', '{0}.0000.h5'.format(base)), 'r')
        chkpt1 = h5py.File (os.path.join('data', '{0}.0001.h5'.format(base)), 'r')
        d0 = chkpt0['primitive']['density'][...]
        d1 = chkpt1['primitive']['density'][...]

        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax1.plot(d0, '-o', mfc='none', label=r'$\rho(t=0)$')
        ax1.plot(d1, '-o', mfc='none', label=r'$\rho(t=1)$')
        ax1.legend(loc='best')
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

    plotters = [Shocktube1(), Shocktube3(), DensityWave()]

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