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
dir_n = 'Drosophila_mCherry_Membrane_02/beads/2016-01-20_16.14.22/Stack_0_Channel_0/'
dir_m = 'Drosophila_Membrane02/Stack1/'
mydir = dir_n
names = {'Cam_Left_00000','Cam_Right_00000'}
#names = {'Beads_A000_T00000','Beads_A090_T00000','Beads_A180_T00000','Beads_A270_T00000'}
#names = {'Cam_Right_00031', 'Cam_Right_00032', 'Cam_Right_00033', 'Cam_Right_00034'}
for name in names:
    print name
    file_name =dir + mydir + name + '.h5'
    #print_hdf5_file_structure(file_name)
    file =  h5py.File(file_name, 'r')
    data = file['/Data'][()]

    print 'data.shape=', data.shape
    num_frames = data.shape[0]
    print 'num_of_frames:', num_frames

    path = mydir+name+'/'
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    print ('max =', data.max(),', min' , data.min())
    for i in range(num_frames):
        #im = Image.fromarray(data, 'I;16')
        im = Image.fromarray(data[i][:][:], 'I;16')
        im.save(path+str(i)+'.tiff')
print 'Done'
