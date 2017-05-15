import matplotlib.pyplot as plt
import h5py
import numpy as np


def main():
    chkpt0 = h5py.File("chkpt.0000.h5", "r")
    chkpt1 = h5py.File("chkpt.0001.h5", "r")
    P0 = chkpt0["primitive"][:,0,0,0]
    P1 = chkpt1["primitive"][:,0,0,0]
    plt.plot(P0, 'o', mec='k', mfc='none')
    plt.plot(P1, 's', mec='b', mfc='none')


    xt = chkpt1["t"].value
    x = np.linspace(-1, 1, P0.size)
    y = np.exp (-(x - xt)**2 / 0.025)
    plt.plot(y)

    plt.show()


if __name__ == "__main__":
    main()
