

class MethodOfLinesTVD(object):

    def plot(self):
        figs = [ ]
        for scheme in ['plm2']:
            try:
                figs.append(self.with_scheme(scheme))
            except IOError as e:
                print e, "for scheme", scheme
        return figs

    def with_scheme(self, scheme):
        import os
        import h5py
        import matplotlib.pyplot as plt

        base = '{0}-{1}'.format(self.name, scheme)
        chkpt = h5py.File (os.path.join('data', '{0}.0001.h5'.format(base)), 'r')
        d = chkpt['primitive']['density'][...]
        p = chkpt['primitive']['pressure'][...]
        b1 = chkpt['primitive']['magnetic1'][...]
        b2 = chkpt['primitive']['magnetic2'][...]
        b3 = chkpt['primitive']['magnetic3'][...]
        m = chkpt['diagnostic']['monopole'][...]
        q = (b1**2 + b2**2 + b3**2) / 2

        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        cmd = ax1.imshow(d, interpolation='nearest')
        cmp = ax2.imshow(p, interpolation='nearest')
        cmq = ax3.imshow(q, interpolation='nearest')
        cmm = ax4.imshow(m, interpolation='nearest')
        ax1.set_title(r'$\rho$')
        ax2.set_title(r'$p$')
        ax3.set_title(r'$B^2$')
        ax4.set_title(r'$\nabla \cdot \vec B$')
        fig.colorbar(cmd, ax=ax1)
        fig.colorbar(cmp, ax=ax2)
        fig.colorbar(cmq, ax=ax3)
        fig.colorbar(cmm, ax=ax4)

        for ax in [ax1, ax2, ax3, ax4]:
	        ax.set_xticks([])
	        ax.set_yticks([])

        fig.suptitle(base)
        return fig


class CylindricalBlast(MethodOfLinesTVD): name = "CylindricalBlast"
class FieldLoop(MethodOfLinesTVD): name = "FieldLoop"
class AbcField(MethodOfLinesTVD): name = "AbcField"


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

    plotters = [CylindricalBlast(), FieldLoop(), AbcField()]

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
