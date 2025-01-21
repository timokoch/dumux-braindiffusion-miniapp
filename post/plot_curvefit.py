# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt
import numpy as np
import glob

# Data for plotting
t = np.genfromtxt("times.dat")
pt = np.genfromtxt("plottimes.dat")

# glob all files with the pattern curvefit-*.dat
files = sorted(glob.glob("curvefit-*.dat"))
files_data = sorted(glob.glob("data-*.dat"))

batch_size = 100
for i in range(0, len(files), batch_size):
    fig, axs = plt.subplots(10, 10, figsize=(15, 9), sharex=True, sharey=True)
    for j, (f, d) in enumerate(zip(files[i:i+batch_size], files_data[i:i+batch_size])):
        data = np.genfromtxt(f)
        data_data = np.genfromtxt(d)
        ax = axs[j // 10, j % 10]
        ax.plot(t/86400.0, data_data, 'ro', label='data')
        ax.plot(pt/86400.0, data, 'b-', label='fit')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
