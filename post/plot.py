#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import click

sec_to_days = 1.0/(24*60*60)

# plot time series data from total.txt, whitematter.txt, and graymatter.txt
@click.command()
@click.option('--folder', default="../build-cmake/app", help='Folder containing the data')
@click.option('-c', '--plot_concentration', is_flag=True, default=False, help='Wheather to plot concentration instead of amount')
def plot_data(folder: str, plot_concentration: bool):
    # load data
    if plot_concentration:
        total = np.loadtxt(f'{folder}/total.txt')
        white = np.loadtxt(f'{folder}/whitematter.txt')
        gray = np.loadtxt(f'{folder}/graymatter.txt')

        total_data = np.loadtxt(f'{folder}/total_data.txt')
        white_data = np.loadtxt(f'{folder}/whitematter_data.txt')
        gray_data = np.loadtxt(f'{folder}/graymatter_data.txt')
    else:
        total = np.loadtxt(f'{folder}/total_amount.txt')
        white = np.loadtxt(f'{folder}/whitematter_amount.txt')
        gray = np.loadtxt(f'{folder}/graymatter_amount.txt')

        total_data = np.loadtxt(f'{folder}//total_amount_data.txt')
        white_data = np.loadtxt(f'{folder}/whitematter_amount_data.txt')
        gray_data = np.loadtxt(f'{folder}/graymatter_amount_data.txt')

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
