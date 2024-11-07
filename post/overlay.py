import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import click
from pathlib import Path


def overlay_images(overlay_data: np.ndarray, background_data: np.ndarray, ax, slicer_op, min_val=None, max_val=None):
    # Choose a specific slice index (for example, slice along the z-axis)
    background_slice = slicer_op(background_data)
    overlay_slice = slicer_op(overlay_data)

    # Create a mask where overlay_slice has NaN values
    mask = np.isnan(overlay_slice)

    # Plot the background image in grayscale
    ax.imshow(background_slice, cmap='gray', interpolation='none')

    # Overlay the second image where it's not NaN, using a colormap and some transparency
    # Masked array so that NaNs in overlay are transparent
    overlay_masked = np.ma.masked_where(mask, overlay_slice)
    ax.imshow(overlay_masked, cmap='inferno', interpolation='none', alpha=1.0, vmin=min_val, vmax=max_val)

    # Display the final overlay image
    ax.axis("off")


def niifile_to_ndarray(nii: Path) -> np.ndarray:
    return nib.load(nii).get_fdata()


def overlay_image_files(overlay: Path, background: Path, ax, slicer_op, min_val=None, max_val=None):
    # Load the two NIfTI images and extract the image data
    overlay_images(
        background_data=niifile_to_ndarray(background),
        overlay_data=niifile_to_ndarray(overlay),
        ax=ax, slicer_op=slicer_op,
        min_val=min_val, max_val=max_val,
    )


@click.command()
@click.option("--overlay", type=Path, required=True)
@click.option("--background", type=Path, required=True)
def plot_data(overlay: Path, background: Path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    overlay_image_files(overlay, background, ax)
    plt.show()

if __name__ == "__main__":
    plot_data()