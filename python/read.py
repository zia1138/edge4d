import h5py
import numpy as np
from PIL import Image
import sys
import os
import errno
import glob


def print_hdf5_file_structure(file_name):
    file = h5py.File(file_name, 'r')
    item = file #["/Configure:0000/Run:0000"]
    print_hdf5_item_structure(item)
    file.close()

def print_hdf5_item_structure(g, offset = '    '):
    if isinstance(g, h5py.File):
        print (g.file, '(File)', g.name)

    elif isinstance(g, h5py.Dataset):
        print ('(Dataset)', g.name, '    len =', g.shape)

    elif isinstance(g, h5py.Group):
        print ('(Group)', g.name)

    else:
        print ('Warning: unknown item in file', g.name)
        sys.exit("EXECUTION IS TERMINATED")

    if isinstance(g, h5py.File) or isinstance(g, h5py.Group):
        for key, val in dict(g).items():
            subg = val
            print (offset, key)
            print_hdf5_item_structure(subg, offset+'    ')


def extract_and_save(file_name, dst_dir):
    try:
        os.makedirs(dst_dir)
    except OSError:
        if not os.path.isdir(dst_dir):
            raise
    
    print(file_name)
    print_hdf5_file_structure(file_name)
    file =  h5py.File(file_name, 'r')
    data = file['/Data'][()]
    thumb = file['/Thumb'][()]

    print ('data.shape=', data.shape)
    num_frames = data.shape[0]
    print ('num_of_frames:', num_frames)

    for i in range(num_frames):
        im = Image.fromarray(data[i][:][:], 'I;16')
        outfilename = dst_dir + "/" + ("%5.5d" % i ) + ".tiff"
        print(outfilename)
        im.save(outfilename)


        
#srcdir = '../../SPIMDataTest/Drosophila_Membrane02/'

#srcdir = "../../SPIMDataTest/Drosophila_mCherry_Membrane_02/2016-01-19_17.47.58/"
#srcdir = "../../SPIMDataTest/Drosophila_mCherry_Membrane_02/beads/2016-01-20_16.09.38/"
srcdir = "../../SPIMDataTest/Drosophila_mCherry_Membrane_02/beads/2016-01-20_16.14.22/"

stack01dirs = ["Stack_0_Channel_0", "Stack_1_Channel_0"]

for subdir in stack01dirs:
    camfiles = glob.glob(srcdir + subdir +  "/*_00000.h5")
    for cam in camfiles:
        camname = cam.split("/")[-1:]
        dst_dir = subdir + "/" + camname[0]
        extract_and_save(cam, dst_dir)

    
for subdir in stack01dirs:
    camfiles = glob.glob(subdir +  "/*_00000.h5")
    for cam in camfiles:
        camname = cam.split("/")[-1:][0]
        files =glob.glob(cam + "/*.tiff")
        files = sorted(files)
        dst_tiff = subdir + "_" + camname.replace("h5", "tiff")
        files = " ".join(files)
        os.system("tiffcp " + files + " " + dst_tiff)

#final_tiff = glob.glob("Stack_[0,1]_Channel_[0,1]/*.tiff")
#final_tiff = sorted(final_tiff)

#os.system("tiffcp " + " ".join(final_tiff) + " final.tiff")
