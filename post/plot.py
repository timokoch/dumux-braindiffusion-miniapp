# SPDX-FileCopyrightText: Copyright © Jørgen Riseth and Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

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
        cgray = np.loadtxt(f'{folder}/corticalgraymatter.txt')
        scgray = np.loadtxt(f'{folder}/subcorticalgraymatter.txt')

        total_data = np.loadtxt(f'{folder}/total_data.txt')
        white_data = np.loadtxt(f'{folder}/whitematter_data.txt')
        cgray_data = np.loadtxt(f'{folder}/corticalgraymatter_data.txt')
        scgray_data = np.loadtxt(f'{folder}/subcorticalgraymatter_data.txt')
    else:
        total = np.loadtxt(f'{folder}/total_amount.txt')
        white = np.loadtxt(f'{folder}/whitematter_amount.txt')
        cgray = np.loadtxt(f'{folder}/corticalgraymatter_amount.txt')
        scgray = np.loadtxt(f'{folder}/subcorticalgraymatter_amount.txt')

        total_data = np.loadtxt(f'{folder}//total_amount_data.txt')
        white_data = np.loadtxt(f'{folder}/whitematter_amount_data.txt')
        cgray_data = np.loadtxt(f'{folder}/corticalgraymatter_amount_data.txt')
        scgray_data = np.loadtxt(f'{folder}/subcorticalgraymatter_amount_data.txt')

    # plot simulation data
    pt = plt.plot(total[:, 0]*sec_to_days, total[:, 1], label='Total')
    pw = plt.plot(white[:, 0]*sec_to_days, white[:, 1], label='White Matter')
    pg = plt.plot(cgray[:, 0]*sec_to_days, cgray[:, 1], label='Cortical Gray Matter')
    pg2 = plt.plot(scgray[:, 0]*sec_to_days, scgray[:, 1], label='Subcortical Gray Matter')

    # plot MRI data
    plt.plot(total_data[:, 0]*sec_to_days, total_data[:, 1], 'x', color=pt[0].get_color(), label='Total Data')
    plt.plot(white_data[:, 0]*sec_to_days, white_data[:, 1], 'x', color=pw[0].get_color(), label='White Matter Data')
    plt.plot(cgray_data[:, 0]*sec_to_days, cgray_data[:, 1], 'x', color=pg[0].get_color(), label='Cortical Gray Matter Data')
    plt.plot(scgray_data[:, 0]*sec_to_days, scgray_data[:, 1], 'x', color=pg2[0].get_color(), label='Subcortical Gray Matter Data')

    plt.xlabel('Time (days)')
    if plot_concentration:
        plt.ylabel('Concentration (mmol/L)')
    else:
        plt.ylabel('Amount (mmol)')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_data()
