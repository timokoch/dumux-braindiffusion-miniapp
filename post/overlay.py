import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import click
from pathlib import Path

def overlay_images(overlay: Path, background: Path, ax, slicer_op):
    # Load your two NIfTI images
    background_nii = nib.load(background)
    overlay_nii = nib.load(overlay)

    # Extract the image data
    background_data = background_nii.get_fdata()
    overlay_data = overlay_nii.get_fdata()

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
    ax.imshow(overlay_masked, cmap='viridis', interpolation='none', alpha=1.0)

    # Display the final overlay image
    ax.axis("off")


@click.command()
@click.option("--overlay", type=Path, required=True)
@click.option("--background", type=Path, required=True)
def plot_data(overlay: Path, background: Path):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    overlay_images(overlay, background, ax)
    plt.show()

if __name__ == "__main__":
    plot_data()