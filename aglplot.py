#!/usr/bin/env python3.5
from pylab import * 
import numpy as np 
import os, sys
import pretty as pr 
import scipy.optimize as opt

'''
aglplot.py

aglplot.py contains a number of plotting routines for deadling
with outputs from the MH3D code.
'''

DPI = 150

cm_use = pr.get_viridis()

def arrow_plot(folder, plane, time, twod=False):

	clf()
	B = np.loadtxt("{}/agl{}/{}_B.agl".format(folder, time, plane))

	dim1 = plane[0]
	dim2 = plane[1]
	bx = np.loadtxt("{}/agl{}/{}_b{}.agl".format(folder, time, plane, dim1))
	by = np.loadtxt("{}/agl{}/{}_b{}.agl".format(folder, time, plane, dim2))

	rho = np.loadtxt("{}/agl{}/{}_rho.agl".format(folder, time, plane))

	pcolormesh(B.T, cmap = cm_use)
	colorbar()

	print (bx.shape)
	number = 0
	for ix in range(bx.shape[0]):
		for iy in range(bx.shape[1]):
			bbx = bx[ix,iy]
			bby = by[ix,iy]

			if twod:
				bbx = by[ix,iy]
				bby = bx[ix,iy]

			ratio = 3 

			dx = ratio * np.log10(np.fabs(bbx)) * np.sign(bbx)
			dy = ratio * np.log10(np.fabs(bby)) * np.sign(bby)

			 

			if B[ix,iy] > (np.mean(B) + (np.std(B)*2)) and B[ix,iy]>1:

				arrow(ix, iy, dy, dx, fc='w', ec='w',head_length=2, head_width=2, alpha=0.5)
				#scatter(ix,iy, alpha=0.5, c="r")
				number += 1 

	print (np.sum( (bbx<0) ))
 	#arrow(0,0,50,50,head_length=2, head_width=2, lw=0.5, fc='w', ec='w')
	xlabel("$x$", fontsize=20)
	ylabel("$y$", fontsize=20)
	xlim(0,256)
	ylim(0,256)
	savefig("arrow_{}_{}.png".format(folder, time), dpi=200)

def correlation_coefficient(patch1, patch2):
    product = np.mean((patch1 - patch1.mean()) * (patch2 - patch2.mean()))
    stds = patch1.std() * patch2.std()
    if stds == 0:
        return 0
    else:
        product /= stds
        return product

def div(run_name, t, plane):

	folder="{}/agl{}".format(run_name, t)
 
	momx = np.loadtxt("{}/{}_momx.agl".format(folder, plane))
	momy = np.loadtxt("{}/{}_momy.agl".format(folder, plane))
	momz = np.loadtxt("{}/{}_momz.agl".format(folder, plane))
	B = np.loadtxt("{}/{}_B.agl".format(folder, plane))
	rho = np.loadtxt("{}/{}_rho.agl".format(folder, plane))
	im = momx + momy + momz

	sh_row, sh_col = im.shape
	correlation = np.zeros_like(im)

	d = 1
	for i in range(d, sh_row - (d + 1)):
		for j in range(d, sh_col - (d + 1)):
				correlation[i, j] = correlation_coefficient(im[i - d: i + d + 1,
                                                        j - d: j + d + 1],
                                                    B[i - d: i + d + 1,
                                                        j - d: j + d + 1])
	#corr = correlation_coefficient(B, im)
	#print (corr)
	#pcolormesh(correlation.T, cmap = cm_use)


	pcolormesh(im.T, cmap = cm_use)
	colorbar()
	contour(B.T, levels=(1,5,10))
	xlim(0,200)
	ylim(0,200)
	#colorbar()

	savefig("corr_{}.png".format(run_name), dpi=100)



def rs10(run_name):

	'''
	make a plot like Riquelme and Spitkovsky 2010
	'''

	# get all the agl folders in the run directory
	folders = os.listdir(run_name)
	rcParams["text.usetex"] = "True"
	nx = 2
	ny = 3

	FONTSIZE = 16

	use_strings=["xz_B", "xz_rho"]
	strings = ["$B$", r"$\rho$"]

	yplot = 0

	times = [2,3,5]

	figure(figsize=(8,8))

	for use_string in use_strings:

		yplot+=1

		for i, tt in enumerate(times):

			folder="agl{}".format(tt)

			folder_path = "{}/{}/".format(run_name, folder)

			print(folder)  

			try:

				files = os.listdir(folder_path)

				for j, f in enumerate(files):

					if use_string in f:
						print (ny,nx,yplot+i*2)
						subplot(ny,nx,yplot+i*2)

						
						position = get_position(ny,nx,yplot+i*2)

						fname = "{}{}".format(folder_path,f)

					 
						d = np.loadtxt(fname )

						print(np.mean(d))


						pcolormesh(d.T, cmap = cm_use)
						colorbar()
						limit = 200
						xlim(0,limit)
						ylim(0,limit)

						if i == 0:
							title(strings[yplot-1], fontsize=20)

						#if "r" in position or mode != "onepage":
						#	cbar.set_label("${}$".format(use_string[3:]), fontsize=FONTSIZE)
						if i == 3:
							xlabel("${}$".format(use_string[0]), fontsize=FONTSIZE)
						else:
							gca().set_xticklabels([])
						if yplot == 1:
							ylabel("${}$".format(use_string[1]), fontsize=FONTSIZE)
						#else:
						#if "l" not in position:
						#	gca().set_yticklabels([])
						#if "b" not in position:
						#	gca().set_xticklabels([])

			except FileNotFoundError:
				print("File not found {}".format(f))

	subplots_adjust(wspace=0.15, hspace =0.1, top=0.96, bottom=0.07)
	savefig("rs10_{}.png".format(run_name), dpi=DPI)
	clf()

	return 0



def file_prep(fname):

	# replace up to first agl string
	# os.system("sed -i -e -n '/agl/,$p' {}".format(fname) )

	# # replace up to first it string
	# os.system("sed -i -e -n '/it/,$p' {}".format(fname) )
	# # substitue first it for i
	# cmd = "awk 'NR==1,/it/{{sub(/it/, \"i\")}} 1' {} > temp.out".format(fname)
	# print (cmd)

	# os.system(cmd)
	# os.system("mv temp.out {}".format(fname) )

	# comment on the lines we don't need
	os.system("sed -i -e 's/agl/# agl/g' {}".format(fname) )
	os.system("sed -i -e 's/it/# it/g' {}".format(fname) )
	os.system("sed -i -e 's/xy/# xy/g' {}".format(fname) )
	os.system("sed -i -e 's/xz/# xz/g' {}".format(fname) )
	os.system("sed -i -e 's/yz/# yz/g' {}".format(fname) )
	os.system("sed -i -e 's/# #/#/g' {}".format(fname) )

	# rm tempfiles
	os.system("rm -f *-e")
	
	return 0

def get_position(ny, nx, i):
	pos = ""
	if ((i-1) % nx) == 0:
		pos += "l"
	elif (i % nx) == 0:
		pos += "r"
	
	if i <= nx:
		pos += "t"
	elif i >= (nx * (ny-1)):
		pos += "b"
	else: pos += "m"

	return pos


def slice_plot(run_name, use_string=None, n=20, mode=None):
	'''
	make a slice plot of a given quantity. 
	The quantity in question is controlled by the use_string.

	Parameters:
		use_string 	string 
					string to plot, e.g. xy_B is the B field
					in an xy slice.
		n 			int
					number of plots to go up to 

	Returns: 0
	'''

	# get all the agl folders in the run directory
	folders = os.listdir(run_name)

	# make a folder to store plots 
	os.system("mkdir slices_{}".format(run_name))

	nx = 3
	ny = (n+2) / nx
	iplot = 1

	if mode == "onepage":
		FONTSIZE = 12
	else:
		FONTSIZE = 20

	if mode == "onepage": 
		#figure(figsize=(3*nx,3*ny))
		suptitle(use_string[0:2]+" "+use_string[3:])

	for i in range(n):

		folder="agl{}".format(i)

		folder_path = "{}/{}/".format(run_name, folder)

		print(folder)  

		try:

			files = os.listdir(folder_path)

			for j, f in enumerate(files):

				if use_string!=None:
					if use_string in f:

						if mode == "onepage":
							subplot(ny,nx,i+1)
						
						position = get_position(ny,nx,i+1)

						fname = "{}{}".format(folder_path,f)

					 
						d = np.loadtxt(fname )

						print(np.mean(d))


						imshow(d.T, cmap = cm_use)
						cbar = colorbar()

						if "r" in position or mode != "onepage":
							cbar.set_label("${}$".format(use_string[3:]), fontsize=FONTSIZE)
						if "b" in position or mode != "onepage":
							xlabel("${}$".format(use_string[0]), fontsize=FONTSIZE)
						if "l" in position or mode != "onepage":
							ylabel("${}$".format(use_string[1]), fontsize=FONTSIZE)

						if "l" not in position and mode == "onepage":
							gca().set_yticklabels([])
						if "b" not in position and mode == "onepage":
							gca().set_xticklabels([])


						if mode != "onepage":
							savefig("slices_{}/t{}{}.png".format(run_name, i,f), dpi=DPI)
							clf()
		except FileNotFoundError:
			print("File not found {}".format(f))

	if mode == "onepage":
		savefig("all{}_{}.png".format(use_string, run_name), dpi=DPI)
		clf()

	return 0

def make_movie(fname="test", use_string = "xy_B", fps = 3):
	'''
	make a movie!
	command is e.g.
	ffmpeg -r 3 -i t%dxy_B.agl.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4
	'''

	cmd = "ffmpeg -r {} -i t%d{}.agl.png -vcodec libx264 -crf 25  -pix_fmt yuv420p {}.mp4".format(fps, use_string, fname)	
	print(cmd)
	isys = os.system(cmd)

	return isys

def fline (x, m, c):
	'''
	a straight line, y = mx+c, for fitting purposes
	'''
	return m*x + c

def db(fnames):

	import astropy.io.ascii as io 
	import pretty

	rcParams["text.usetex"] = "False"

	lstyles = ["-", "--", "-.", "-", "--", "-."]
	colors = pretty.get_colors()

	for i, fname in enumerate(fnames):

		file_prep(fname)

		data = io.read(fname)

		t = data["time"]
		b = (data["bmax"] - data["bmax"][0]) / data["bmax"][0]
		plot(t, b, label="max, {}".format(fname), linewidth=2, ls=lstyles[i], c=colors[0])

		t = data["time"]
		b = (np.sqrt(data["b"]) - np.sqrt(data["b"][0])) / np.sqrt(data["b"][0])
		plot(t, b, label="rms, {}".format(fname), linewidth=2, ls=lstyles[i], c=colors[4])



		pretty.float_legend(loc=2)

	ylim(0.001,1000)
	semilogy()
	xlabel("$t$", fontsize=20)
	ylabel("$\delta B / B_0$", fontsize=20)

	savefig("db.png", dpi=300)

	return 0


def physical(fname1, fname2 = None, fit=True, savename="physical", tt=(2,4)):

	'''
	make a plot of how the physical parameters of the simulation evolves 
	over time
	'''

	import astropy.io.ascii as io 
	import pretty

	pretty.set_pretty()
	rcParams["text.usetex"] = "False"

	# figure out whether we want to plot multiple things
	if fname2 != None:
		fnames = [fname1, fname2]

	else: fnames = [fname1]

	lstyles = ["-", "--"]
	colors = pretty.get_colors()

	for i, fname in enumerate(fnames):

		# run some sed commands on the output file to comment out
		# lines we don't want
		file_prep(fname)

		# read using astropy table format
		data = io.read(fname)

		t = data["time"]

		labels = ["b","ke","thermal"]


		# array to store total energy density
		total = np.zeros(len(data["b"]))

		# now loop over different forms of energy
		icolor = 0
		for l in labels:
			plot(t, data[l], label=l+" "+fname, linewidth=1.5, ls=lstyles[i], c=colors[icolor])
			total += data[l] - data[l][0]
			#total += data[l]
			icolor+=1

		# plot the total and add a legend
		total = total - total[0]
		plot(t, total, label="total - E0", linewidth=3, ls=lstyles[i], c=colors[icolor])
		pretty.float_legend(loc=4)


		if fit:
			# now fit the ln(total) line 
			x = t 
			y = np.log(total)
			t1 = np.argmin(np.fabs(t - tt[0]))
			t2 = np.argmin(np.fabs(t - tt[1]))
			print (len(x[t1:t2]), t1, t2, tt, t)
			popt1 = opt.curve_fit(fline, x[t1:t2], y[t1:t2])[0] 
			print(popt1)
			
			# plot the fitted line by converting back to ln space
			plot(x, np.exp(fline(x, popt1[0], popt1[1])), c="k", linestyle="--", label=None)

			# repeat process in linear regime
			#t1 = np.argmin(np.fabs(t - 10))
			#t2 = np.argmin(np.fabs(t - 15))
			#popt2 = opt.curve_fit(fline, x[t1:t2], y[t1:t2])[0] 
			#print(popt2)
			#plot(x, np.exp(fline(x, popt2[0], popt2[1])), c="k", linestyle="--", label=None)

			#print(popt2/popt1)
			#print(popt2/2.512)

			text(2,100, "Rate: {:.2f}".format(popt1[0]))


	ylim(0.1,1000)
	xlim(0,22)
	semilogy()
	#pretty.float_legend(loc=4)
	xlabel("t", fontsize=20)
	ylabel("E", fontsize=20)

	savefig("growth_{}.png".format(savename), dpi=DPI)

	return 0


if __name__ == "__main__":

	mode = sys.argv[1]

	if mode == "slice":

		for j in range(2, len(sys.argv[:-1])):

			run_name = sys.argv[j]

			# make the slice plots
			use_string=sys.argv[-1]

			if use_string == "all":
				slice_plot(run_name, use_string=None)
			else:
				slice_plot(run_name, use_string=use_string)

	elif mode == "one":

		for j in range(2, len(sys.argv[:-2])):

			run_name = sys.argv[j]

			# make the slice plots
			use_string=sys.argv[-2]
			n = int(sys.argv[-1])

			if use_string == "all":
				slice_plot(run_name, use_string=None, mode="onepage", n=n)
			else:
				slice_plot(run_name, use_string=use_string, mode="onepage", n=n)

	elif mode == "mov":

		use_string=sys.argv[-1]

		# combine plots into movie with 3 frames per second
		make_movie(fname = "test", use_string=use_string, fps=3)

	elif mode == "t":

		run_names = sys.argv[2:]
		over_time(run_names)

	elif mode == "p":

		fname = sys.argv[2]
		savename = sys.argv[3]
		times = (float(sys.argv[4]), float(sys.argv[5]))

		physical(fname, savename=savename, tt=times)

	elif mode == "db":
		run_names = sys.argv[2:]
		db(run_names)

	elif mode == "pn":

		fname = sys.argv[2]
		savename = sys.argv[3]

		physical(fname, fit=False, savename=savename)

	elif mode == "div":
		fname = sys.argv[2]
		t = int(sys.argv[3])
		plane = sys.argv[4]

		div(fname, t, plane)

	elif mode == "pcomp":

		fname = sys.argv[2]
		fname2 = sys.argv[3]

		physical(fname, fname2, savename="comp", fit=False)

	elif mode == "rs10":

		run_name = sys.argv[2]

		rs10(run_name)

	elif mode == "arrow":
		folder = sys.argv[2]
		time = int(sys.argv[3])
		plane = sys.argv[4]

		for t in range(time):
			arrow_plot(folder, plane, t)

	elif mode == "a2d":
		folder = sys.argv[2]
		time = int(sys.argv[3])
		plane = sys.argv[4]

		for t in range(time):
			arrow_plot(folder, plane, t, twod=True)



