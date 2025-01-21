# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Script that maps the vtk simulation data to the MRI image space
"""

from vtk2mri import pyvista_mesh_to_mri
import xml.etree.ElementTree as ET
import multiprocessing as mp
import pyvista as pv
import numpy as np
import simple_mri as sm

# file prefix (hard-coded for now)
PREFIX = "../build-cmake/app/"

# time stamps at which we have data
TIME_STAMPS = [0, 15716, 86001, 172280, 251641]
VTK_FIELDS = ["c", "c_data"]

# extract data from the VTK files and convert to MRI
# this is a CPU-bound task, so we use multiprocessing to speed up the process
def process_vtk_file(vtk_file, ref_mri, time_stamp):
    mesh = pv.read(f'{PREFIX}{vtk_file}')
    newmris = pyvista_mesh_to_mri(mesh, VTK_FIELDS, ref_mri)
    for data_name in VTK_FIELDS:
        file_name = f'{PREFIX}{data_name}_{time_stamp:d}.nii.gz'
        print(f'Saving mapped field {data_name} to {file_name}')
        sm.save_mri(newmris[data_name], file_name, np.single)

def get_vtk_files(root):
    # get the names of files (attribute "file") with timesteps (attribute "timestep") from the PVD file
    # each file is a DataSet tag in the PVD file
    vtk_files = []
    for data_set in root.findall('.//DataSet'):
        timestep = int(data_set.attrib['timestep'])
        if timestep not in TIME_STAMPS:
            continue
        vtk_file = data_set.attrib['file']
        print(f'Found simulation output file {vtk_file} @t={timestep}')
        vtk_files.append(vtk_file)
    return vtk_files

if __name__ == '__main__':
    # read the PVD file which is an XML file that contains the list of VTK files
    vtk_files = get_vtk_files(ET.parse(f'{PREFIX}diffusion.pvd').getroot())
    ref_filename = '../data/sub-01_ses-01_T1w_registered.nii.gz'
    ref_mri = sm.load_mri(ref_filename, dtype=np.single)
    print(f'Loaded reference MRI image from {ref_filename}')

    print("Mapping simulation data to MRI... (this may take a while)")
    with mp.Pool(len(TIME_STAMPS)) as pool:
        pool.starmap(process_vtk_file, [(vtk_file, ref_mri, TIME_STAMPS[i]) for i, vtk_file in enumerate(vtk_files)])
