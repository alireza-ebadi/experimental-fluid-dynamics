import numpy as np
import pandas as pd
from scipy.stats import linregress
import glob
import tkinter as tk            # for the gui
from tkinter import filedialog
from numDeri import numericalDerivative as deri

#############################################
### Function sorts and processes PIV data and outputs the velocity profiles
def PIVpro(xoffset, yoffset, xcoor, ycoor):

	#############################################
	### Function sorts and processes PIV data
	## inputs:
	# xoffset (mm): The x offset from the origin
    # yoffset (mm): The wall location in images
    # xcoor (-): = 1 if the flow and the x-coordinate have the same sign, i.e flow
    # goes from left to right; otherwise -1
    # ycoor (-): = 1 if y direction is from bottom to top; otherwise -1
    ## outputs:
    # U (m/s): ensemble average streamwise velocity: 2D matrix (columns = streamwise position, rows = wall-normal position)
    # V (m/s): ensemble average wall-normal velocity: 2D matrix
    # urms (m/s): rms streamwise velocity: 2D matrix
    # vrms (m/s): rms wall-normal velocity: 2D matrix
    # uv (m2/s2): Reynolds stress in wall-normal direction: 2D matrix
    # y (m): wall-normal grid points
    # x (m): streamwise grid points
	# uMaster (m/s): master streamwise velocity matrix: 3D matrix (columns = streamwise position, rows = wall-normal position, 3rd dimension = snapshots)
    # vMaster (m/s): master wall-normal velocity matrix: 3D matrix
	# velMaster (m/s): master velocity matrix: 2D matrix (columns = snapshots, first half of the rows = streamwise velocity, 2nd half of the rows = wall-normal velocity)
	# duMasterdy (1/s): duMaster/dy
	#############################################

	# for the gui window
	root = tk.Tk()
	root.withdraw()
	root.attributes("-topmost", True)

	# select the main directory
	Maindir = filedialog.askdirectory()

    # count number of text files in the directory (number of snapshots)
	nSnaps = len(glob.glob1(Maindir,"*.txt"))

	# read a sample text file to find the grid size
	temp = pd.read_table(Maindir + '/B00001.txt',  sep='\s+')
	temp = np.array(temp)
	temp = temp[:, :4]

    # number of grid points in x, and y
	nGridx = len(np.unique(temp[:, 0]))
	nGridy = len(np.unique(temp[:, 1]))
	nGridxy = nGridx*nGridy

    # master matrix for x, and y converted to meters
	xMaster = (temp[:, 0] - xoffset)*1e-3
	yMaster = (temp[:, 1] - yoffset)*1e-3

    # unique arrays fro x and y
	x = xMaster[:nGridx]
	y = yMaster[:nGridxy:nGridx]

    # flag negative values
	x[x<0] = np.nan
	y[y<0] = np.nan

    # create master velocity matrices
	uMaster = np.zeros((nGridy, nGridx, nSnaps))
	vMaster = np.zeros((nGridy, nGridx, nSnaps))
	velMaster = np.zeros((nGridxy*2, nSnaps))
	duMasterdy = np.zeros((nGridy, nGridx, nSnaps))


	# velocity below threshold is considered 0
	thresh = 1e-8
	for s in range(nSnaps):
		temp = pd.read_table(Maindir + '/B' + str(s + 1).zfill(5) + '.txt',  sep='\s+')
		temp = np.array(temp)
		u = xcoor*temp[:, 2]
		u[np.abs(u) < thresh] = np.nan
		v = ycoor*temp[:, 3]
		v[np.abs(v) < thresh] = np.nan
		uMaster[:, :, s] = np.reshape(u, (nGridy, nGridx))
		vMaster[:, :, s] = np.reshape(v, (nGridy, nGridx))
		velMaster[:, s] = np.concatenate((u, v), axis=0)
		for i in range(nGridx):
			utemp = uMaster[:, i, s]
			dudy = deri(y, utemp)
			duMasterdy[:, i, s] = dudy


    # flip matrices so they read from bottom to top
	y = y[::-1]
	uMaster = uMaster[::-1, :, :]
	vMaster = vMaster[::-1, :, :]
	duMasterdy = duMasterdy[::-1, :, :]


    # statistics
	U = np.nanmean(uMaster, 2)
	V = np.nanmean(vMaster, 2)
	uprime = np.zeros((nGridy, nGridx, nSnaps))
	vprime = np.zeros((nGridy, nGridx, nSnaps))
	for s in range(nSnaps):
		uprime[:,:, s] = uMaster[:,:, s] - U
		vprime[:,:, s] = vMaster[:,:, s] - V
	urms = np.nanmean(uprime**2, 2)**0.5
	vrms = np.nanmean(vprime**2, 2)**0.5
	uvprime = uprime*vprime
	uv = np.nanmean(uvprime, 2)
	return U, V, urms, vrms, uv, y, x, velMaster, uMaster, vMaster, duMasterdy

def PIVpro2(xoffset, yoffset, xcoor, ycoor, nstd, thresh, Maindir):

	#############################################
	### Function sorts and processes PIV data
	## inputs:
	# xoffset (mm): The x offset from the origin
    # yoffset (mm): The wall location in images
    # xcoor (-): = 1 if the flow and the x-coordinate have the same sign, i.e flow
    # goes from left to right; otherwise -1
    # ycoor (-): = 1 if y direction is from bottom to top; otherwise -1
	# nstad (-): filter elements larger than mean(u) +/- nstd*std(u)
	# thresh (m/s): elements smaller than thresh are considered 0
    ## outputs:
    # U (m/s): ensemble average streamwise velocity: 2D matrix (columns = streamwise position, rows = wall-normal position)
    # V (m/s): ensemble average wall-normal velocity: 2D matrix
    # urms (m/s): rms streamwise velocity: 2D matrix
    # vrms (m/s): rms wall-normal velocity: 2D matrix
    # uv (m2/s2): Reynolds stress in wall-normal direction: 2D matrix
    # y (m): wall-normal grid points
    # x (m): streamwise grid points
	# uMaster (m/s): master streamwise velocity matrix: 3D matrix (columns = streamwise position, rows = wall-normal position, 3rd dimension = snapshots)
    # vMaster (m/s): master wall-normal velocity matrix: 3D matrix
	# velMaster (m/s): master velocity matrix: 2D matrix (columns = snapshots, first half of the rows = streamwise velocity, 2nd half of the rows = wall-normal velocity)
	# duMasterdy (1/s): duMaster/dy
	#############################################

	# for the gui window
	#root = tk.Tk()
	#root.withdraw()
	#root.attributes("-topmost", True)

	# select the main directory
	#Maindir = filedialog.askdirectory()

    # count number of text files in the directory (number of snapshots)
	nSnaps = len(glob.glob1(Maindir,"*.txt"))

	# read a sample text file to find the grid size
	temp = pd.read_table(Maindir + '/B00001.txt',  sep='\s+')
	temp = np.array(temp)
	temp = temp[:, :4]

    # number of grid points in x, and y
	nGridx = len(np.unique(temp[:, 0]))
	nGridy = len(np.unique(temp[:, 1]))
	nGridxy = nGridx*nGridy

    # master matrix for x, and y converted to meters
	xMaster = (temp[:, 0] - xoffset)*1e-3
	yMaster = (temp[:, 1] - yoffset)*1e-3

    # unique arrays fro x and y
	x = xMaster[:nGridx]
	y = yMaster[:nGridxy:nGridx]

    # flag negative values
	x[x<0] = np.nan
	y[y<0] = np.nan

    # create master velocity matrices
	uMaster = np.zeros((nGridy, nGridx, nSnaps))
	vMaster = np.zeros((nGridy, nGridx, nSnaps))
	velMaster = np.zeros((nGridxy*2, nSnaps))
	duMasterdy = np.zeros((nGridy, nGridx, nSnaps))


	# velocity below threshold is flagged
	for s in range(nSnaps):
		temp = pd.read_table(Maindir + '/B' + str(s + 1).zfill(5) + '.txt',  sep='\s+')
		temp = np.array(temp)
		u = xcoor*temp[:, 2]
		u[np.abs(u) < thresh] = np.nan
		v = ycoor*temp[:, 3]
		v[np.abs(v) < thresh] = np.nan
		uMaster[:, :, s] = np.reshape(u, (nGridy, nGridx))
		vMaster[:, :, s] = np.reshape(v, (nGridy, nGridx))
		velMaster[:, s] = np.concatenate((u, v), axis=0)
		for i in range(nGridx):
			utemp = uMaster[:, i, s]
			dudy = deri(y, utemp)
			duMasterdy[:, i, s] = dudy


    # flip matrices so they read from bottom to top
	y = y[::-1]
	uMaster = uMaster[::-1, :, :]
	vMaster = vMaster[::-1, :, :]
	duMasterdy = duMasterdy[::-1, :, :]


    # statistics
	U = np.nanmean(uMaster, 2)
	V = np.nanmean(vMaster, 2)

	# standard deviation filter: remove elements larger than mean(u) + nstd*std(u)
	ustd = np.nanstd(uMaster, 2)
	for s in range(nSnaps):
		utemp = uMaster[:, :, s]
		utemp[np.abs(utemp)>np.abs(U + nstd*ustd)] = np.nan;
		uMaster[:, :, s] = utemp

		vtemp = vMaster[:, :, s]
		vtemp[np.abs(utemp)>np.abs(U + nstd*ustd)] = np.nan;
		vMaster[:, :, s] = vtemp


	uprime = np.zeros((nGridy, nGridx, nSnaps))
	vprime = np.zeros((nGridy, nGridx, nSnaps))
	for s in range(nSnaps):
		uprime[:,:, s] = uMaster[:,:, s] - U
		vprime[:,:, s] = vMaster[:,:, s] - V
	urms = np.nanmean(uprime**2, 2)**0.5
	vrms = np.nanmean(vprime**2, 2)**0.5
	uvprime = uprime*vprime
	uv = np.nanmean(uvprime, 2)
	return U, V, urms, vrms, uv, y, x, velMaster, uMaster, vMaster, duMasterdy
