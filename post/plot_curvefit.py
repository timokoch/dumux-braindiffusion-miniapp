# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt
import numpy as np
import glob
import click

@click.command()
@click.option('--folder', default="../build-cmake/app", help='Folder containing the data')
@click.option('-n', "--number-of-plots", default=2, help='Number of plots before stopping (integer)', type=int)
def plot_curves(folder: str, number_of_plots: int):
    t = np.genfromtxt(f"{folder}/times.dat")
    pt = np.genfromtxt(f"{folder}/plottimes.dat")

    # glob all files with the pattern curvefit-*.dat
    files = sorted(glob.glob(f"{folder}/curvefit-*.dat"))
    files_data = sorted(glob.glob(f"{folder}/data-*.dat"))

    batch_size = 100
    for i in range(0, len(files), batch_size):
        if i >= number_of_plots*batch_size:
            break

        fig, axs = plt.subplots(10, 10, figsize=(15, 9), sharex=True, sharey=True)
        for j, (f, d) in enumerate(zip(files[i:i+batch_size], files_data[i:i+batch_size])):
            data = np.genfromtxt(f)
            data_data = np.genfromtxt(d)
            ax = axs[j // 10, j % 10]
            ax.plot(t/86400.0, data_data, 'ro', label='data')
            ax.plot(pt/86400.0, data, 'b-', label='fit')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0, hspace=0)
        plt.show()


if __name__ == '__main__':
    plot_curves()
