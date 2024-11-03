#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

sec_to_days = 1.0/(24*60*60)

plot_concentration = True

# plot time series data from total.txt, whitematter.txt, and graymatter.txt
def plot_data():
    # load data
    if plot_concentration:
        total = np.loadtxt('total.txt')
        white = np.loadtxt('whitematter.txt')
        gray = np.loadtxt('graymatter.txt')

        total_data = np.loadtxt('total_data.txt')
        white_data = np.loadtxt('whitematter_data.txt')
        gray_data = np.loadtxt('graymatter_data.txt')
    else:
        total = np.loadtxt('total_amount.txt')
        white = np.loadtxt('whitematter_amount.txt')
        gray = np.loadtxt('graymatter_amount.txt')

        total_data = np.loadtxt('total_amount_data.txt')
        white_data = np.loadtxt('whitematter_amount_data.txt')
        gray_data = np.loadtxt('graymatter_amount_data.txt')

    # plot simulation data
    pt = plt.plot(total[:, 0]*sec_to_days, total[:, 1], label='Total')
    pw = plt.plot(white[:, 0]*sec_to_days, white[:, 1], label='White Matter')
    pg = plt.plot(gray[:, 0]*sec_to_days, gray[:, 1], label='Gray Matter')

    # plot MRI data
    plt.plot(total_data[:, 0]*sec_to_days, total_data[:, 1], 'x', color=pt[0].get_color(), label='Total Data')
    plt.plot(white_data[:, 0]*sec_to_days, white_data[:, 1], 'x', color=pw[0].get_color(), label='White Matter Data')
    plt.plot(gray_data[:, 0]*sec_to_days, gray_data[:, 1], 'x', color=pg[0].get_color(), label='Gray Matter Data')

    plt.xlabel('Time (days)')
    if plot_concentration:
        plt.ylabel('Concentration (mmol/L)')
    else:
        plt.ylabel('Amount (mmol)')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_data()
