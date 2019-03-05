from __future__ import division

import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import viz
import time
import copy
import data_utils
import scipy.io

def main():
    def read_mat(mat_path):
        mat=scipy.io.loadmat(mat_path)
        #a=mat['position'].transpose(2,1,0)
        a=mat['position']
        return a # ndarray type

    xyz_info = read_mat(file) # xyz,joint,frame
    frame_number = xyz_info.shape(2) # get frame number of current file
    xyz_gt = xyz_info.transpose(2,1,0)


    # === Plot and animate ===

    fig = plt.figure()
    ax = plt.gca(projection='3d')
    ob = viz.Ax3DPose(ax)

    # Plot the conditioning ground truth
    for i in range(frame_number):
        ob.update(xyz_gt[i, :])
        plt.show(block=False)
        fig.canvas.draw()
        plt.pause(0.01)


if __name__ == '__main__':
    main()
