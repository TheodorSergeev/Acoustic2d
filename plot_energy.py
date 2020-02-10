import os
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import csv

#https://www.epcc.ed.ac.uk/blog/2016/08/23/mpi-performance-knl

def plot_execution_time(filename):
    energy = []
    time   = []
    
    with open(filename, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for row in plots:
            energy.append(float(row[0]))
            time  .append(float(row[1]))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.27, 11.69))

    tickssize = 10
    labelsize = 12
    titlesize = 14
    suptlsize = 18

    fig.suptitle("", fontsize=suptlsize, fontweight='bold')

    color_x  = "#B72D00"
    color_y  = "#B72FFF"
    color_xy = "#6CDB39"

    params = dict(  lw              = 1.6, 
                    #color           = color,
                    alpha           = 1.0,
                    markersize      = 4,
                    marker          = 'o',
                    markevery       = -1,
                    markerfacecolor = 'white',
                    #markeredgecolor = color,
                    markeredgewidth = 1.2) 

    plt.subplots_adjust(top=0.9, hspace=0.4)

    ax1.tick_params(direction='inout', length=6, width=1.5, colors='black')

    ax1.grid(linestyle=':')

    #ax1.set_xticks(proc_num)
    #ax1.set_xticklabels(msg_size, rotation=40, fontsize=tickssize)

    #ax1.set_ylim(-0.000005, 0.000255)
    ax1.set_xlabel(r't', fontsize=labelsize)
    ax1.set_ylabel(r'Kinetic energy', fontsize=labelsize)
    ax1.set_title (filename, fontsize=titlesize)
    #ax1.plot(proc_num, exec_time_x,  **params, label='x decomposition',  color=color_x)
    ax1.plot(time, energy,  **params, label='asd',  color=color_y)
    #ax1.plot(proc_num, exec_time_xy, **params, label='xy decomposition', color=color_xy)

    ax1.legend()

    plt.savefig(filename + "_plot.pdf", format='pdf')
    #plt.show()


dir_name = "energy"

for filename in os.listdir(dir_name):
    if filename.endswith(".csv"): 
        plot_execution_time(dir_name + "/" + filename)
        continue
    else:
        continue