# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt
import matplotlib as mpl
from overlay import overlay_images, niifile_to_ndarray
from map_to_mri import TIME_STAMPS, VTK_FIELDS, PREFIX
import numpy as np

mpl.rcParams['savefig.pad_inches'] = 0

if __name__ == '__main__':
    fig, axs = plt.subplots(2, len(TIME_STAMPS), figsize=(15, 6), frameon=False)
    background_data = niifile_to_ndarray('../data/mri-dataset-pre-contrast-only/mri_dataset/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii.gz')
    for j in range(len(VTK_FIELDS)):
        for i, time_stamp in enumerate(TIME_STAMPS):
            overlay_images(
                overlay_data=niifile_to_ndarray(f'{PREFIX}{VTK_FIELDS[j]}_{time_stamp}.nii.gz'),
                background_data=background_data,
                ax=axs[j][i],
                #slicer_op=lambda img: img[:, :, 310].T,
                slicer_op=lambda img: np.flipud(img[150, :, :].T),
                min_val=0, max_val=0.05,
            )
            axs[j][i].annotate(
                f'{time_stamp/86400*24:.0f}h',
                xy=(0.04, 0.96), xycoords='axes fraction',
                color='white', fontsize=16,
                ha='left', va='top',
            )

            axs[j][i].margins(0)
            axs[j][i].get_xaxis().set_visible(False)
            axs[j][i].get_yaxis().set_visible(False)

    plt.autoscale(tight=True)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0, left=0.1, right=0.9, bottom=0.1, top=0.9)

    # add vertical colorbar on the right spanning both rows
    cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(axs[0][0].images[1], cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('c in mmol/l', fontsize=18)

    fig.savefig('mapped_mri.png', dpi=600)
    fig.savefig('mapped_mri_small.png', dpi=72)

    plt.show()
