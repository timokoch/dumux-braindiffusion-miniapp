# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: CC0-1.0

import matplotlib.pyplot as plt
import numpy as np

sec_to_days = 1.0/(24*60*60)

# plot time series data from total.txt, whitematter.txt, and graymatter.txt
def plot_data():
    # load data
    total = np.loadtxt('total.txt')
    white = np.loadtxt('whitematter.txt')
    gray = np.loadtxt('graymatter.txt')

    # plot data
    plt.plot(total[:, 0]*sec_to_days, total[:, 1], label='Total')
    plt.plot(white[:, 0]*sec_to_days, white[:, 1], label='White Matter')
    plt.plot(gray[:, 0]*sec_to_days, gray[:, 1], label='Gray Matter')
    plt.xlabel('Time (days)')
    plt.ylabel('Concentration (mmol/L)')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_data()
