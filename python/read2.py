import h5py
import numpy as np
import Image
import sys
import os
import errno

def print_hdf5_file_structure(file_name):
    file = h5py.File(file_name, 'r')
    item = file #["/Configure:0000/Run:0000"]
    print_hdf5_item_structure(item)
    file.close()

def print_hdf5_item_structure(g, offset = '    '):
    if isinstance(g, h5py.File):
        print g.file, '(File)', g.name

    elif isinstance(g, h5py.Dataset):
        print '(Dataset)', g.name, '    len =', g.shape

    elif isinstance(g, h5py.Group):
        print '(Group)', g.name

    else:
        print 'Warning: unknown item in file', g.name
        sys.exit("EXECUTION IS TERMINATED")

    if isinstance(g, h5py.File) or isinstance(g, h5py.Group):
        for key, val in dict(g).iteritems():
            subg = val
            print offset, key
            print_hdf5_item_structure(subg, offset+'    ')
dir = '../../SPIMDataTest/'
dir_n = 'CriWT_H1_TexasRed/'
dir_m = 'Drosophila_Membrane02/'
stacks = ['0', '1']
angles = [0, 90]
mydir = dir_m
for index in [0, 1]:
    stack = stacks[index]
    angle = angles[index]
    #name = 'Cam_Left_00400'
    name = 'Cam_Right_00030'
    file_name =dir + mydir + 'Stack'+stack+'/' + name + '.h5'
    file =  h5py.File(file_name, 'r')
    data = file['/Data'][()]
    #thumb = file['/Thumb'][()]

    print mydir
    print 'data.shape=', data.shape
    num_frames = data.shape[0]
    print 'num_of_frames:', num_frames

    path = mydir+name+'/'
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    for i in range(num_frames):
        #im = Image.fromarray(data, 'I;16')
        im = Image.fromarray(data[i][:][:], 'I;16')
        newname = 'spim_TL'+str(i)+'_Angle'+str(angle)+'.tiff'
        im.save(path+newname)
print 'Done'
