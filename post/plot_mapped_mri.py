import matplotlib.pyplot as plt
from overlay import overlay_images
from map_to_mri import TIME_STAMPS, VTK_FIELDS, PREFIX
import numpy as np

if __name__ == '__main__':
    fig, axs = plt.subplots(2, len(TIME_STAMPS), figsize=(25, 10))
    for j in range(len(VTK_FIELDS)):
        for i, time_stamp in enumerate(TIME_STAMPS):
            overlay_images(
                f'{PREFIX}{VTK_FIELDS[j]}_{time_stamp}.nii.gz',
                '../data/sub-01_ses-01_T1w_registered.nii.gz',
                axs[j][i],
                #lambda img: img[:, :, 310].T,
                lambda img: np.flipud(img[150, :, :].T),
                min_val=0, max_val=0.05,
            )
            axs[j][i].annotate(
                f'{time_stamp/86400*24:.0f}h',
                xy=(0.04, 0.96), xycoords='axes fraction',
                color='white', fontsize=11,
                ha='left', va='top',
            )

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()