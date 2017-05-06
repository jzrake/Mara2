import matplotlib.pyplot as plt
import h5py


def main():
    chkpt0 = h5py.File("chkpt.0000.h5", "r")
    chkpt1 = h5py.File("chkpt.0001.h5", "r")
    P0 = chkpt0["primitive"][...]
    P1 = chkpt1["primitive"][...]
    plt.plot(P0, 'o', mec='k', mfc='none')
    plt.plot(P1, 'o', mec='b', mfc='none')
    plt.show()


if __name__ == "__main__":
    main()
