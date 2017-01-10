import numpy as numpy
from pylab import *
import os, sys
from cobra_sub import smooth
import read_output as rd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

title_string = "mdot1\t\t\t\t\t\t\t\t\t\t\t\tmdot2\t\t\t\t\t\t\t\t\t\t\t\tmdot3"

def create_continuum(cont_points, wave):

	for i in range(len(wave)):
		print( "")

	return flux



def get_linelist(filenames):

	lines = np.array([])

	for i in range(len(filenames)):

		x = np.loadtxt(filenames[i], comments="#", usecols=(1,3), unpack=True)

		lines = np.concatenate((lines, x[1]))

	return lines



def nearest_neighbour(wave, flux):

	last = flux[0]

	for i in range(len(flux)):

		if flux[i] == 0:
			flux[i] = last
		else:
			last = flux[i]

	return flux


def make_f(a, x):

 	deg = len(a)
 	y = 0

 	for i in range(deg):
 		y += a[i] * (x**(deg - i - 1))

 	return y


def get_continuum(wave, spec, lines, lmin=3646, lmax = 6900, deg = 3):

	i_jump = (wave > lmin) * (wave < lmax)

	i_lines = np.zeros(len(wave), dtype=bool)

	for i in range(len(wave)):

		comp = lines - wave[i]
		comp = np.absolute(comp)

		i_line = (comp < 100)

		if np.sum(i_line) == 0:
			i_lines[i] = True
		else:
			i_lines[i] = False		

	i_select = i_jump * i_lines

	continuum_points = wave[i_select]

	#print continuum_points, i_select

	coefficients = np.polyfit(continuum_points, spec[i_select], deg)

	return coefficients



def plot_all (filename):

	#filename = 
	#suffix = "slow2_opt"

	s = rd.read_spec_file(filename)

	rd.setpars()
	plt.rcParams["lines.linewidth"] = 1.5



	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", "data/atomic_macro/he_top_lines.py" ]

	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")



	nspec = len(s.spec)

	figure(figsize=(8.3,8.3))

	for i in range(len(s.spec)):

		subplot(nspec,1,i+1)
		a = get_continuum(s.wavelength, s.spec[i], lines)
		plot(s.wavelength, smooth(s.spec[i]/make_f(a, s.wavelength)) )

		xlim(3100,6900)
		ylim(0.51, np.max(s.spec[i]/make_f(a, s.wavelength)))
		text(6100, 0.7*(np.max(s.spec[i]/make_f(a, s.wavelength)) + 0.51), "$i=%s^{\circ}$" % incs[i])
		if i == 2: ylabel("Flux /continuum")

		ax = gca()

		if i < nspec - 1: ax.set_xticklabels([])

	subplots_adjust(hspace=0)
	xlabel("Wavelength [\AA]")
	savefig("all_%s.png" % filename[:-5])


def make_fig4(filename):

	s = rd.read_spec_file(filename)

	#atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", "data/atomic_macro/he_top_lines.py" ]
	#lines = get_linelist(atomic_filenames)
	#a = get_continuum(s.wavelength, s.spec[i], lines)

	figure(figsize=(8.3,6))

	gs = gridspec.GridSpec(1, 1)
	#subplot(gs[0])
	plot(s.wavelength, 1e12*smooth(s.emitted), label="Total Spectrum")
	#text(6100, 0.5*1e12*np.max(s.emitted), "Total Spectrum")
	ax = gca()
	#ax.set_xticklabels([])
	#ylim(0.01,np.max(s.emitted))

	#subplot(gs[1])
	plot(s.wavelength, 1e12*smooth(s.disk), label="Disk")
	#text(6100, 0.5*1e12*np.max(s.disk), "Disk")
	ax = gca()
	#ax.set_xticklabels([])
	ylabel ("Flux ($10^{-12}$~erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)")
	#ylim(0.01,np.max(s.disk))

	#subplot(gs[2])
	plot(s.wavelength, 1e12*smooth(s.wind), label="Wind")
	#text(6100, 0.5*1e12*np.max(s.wind), "Wind")
	#ylim(0.01,np.max(s.wind))


	xlabel("Wavelength [\AA]")
	subplots_adjust(hspace=0.1)

	savefig('fig4.eps',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()
	#clf()


def make_fig3(filename):

	s = rd.read_spec_file(filename)

	figure(figsize=(8.3,11.6))

	nspec = len(s.spec)

	incs = ("10", "17.5", "45", "62.5", "80")

	for i in range(len(s.spec)):

		subplot(nspec,1,i+1)
		xlim(3100,6900)

		plot(s.wavelength, 1e12*smooth(s.spec[i] ), c="k")

		if i == 2: ylabel("Flux ($10^{-12}$~erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)")

		text(6100, 0.7*1e12*(np.max(s.spec[i])), "$i=%s^{\circ}$" % incs[i])

		ax = gca()
		if i < nspec - 1: ax.set_xticklabels([])

	subplots_adjust(hspace=0)
	xlabel("Wavelength [\AA]")

	savefig('fig3.eps',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()


def make_fig_uv1(filename):

	s = rd.read_spec_file(filename)

	figure(figsize=(8.3,11.6))

	nspec = len(s.spec)

	incs = ("10", "17.5", "45", "62.5", "80")

	for i in range(len(s.spec)):

		subplot(nspec,1,i+1)
		xlim(1000,1800)

		plot(s.wavelength, 1e12*smooth(s.spec[i] ), c="k")

		if i == 2: ylabel("Flux ($10^{-12}$~erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)")

		text(1600, 0.7*1e12*(np.max(s.spec[i])), "$i=%s^{\circ}$" % incs[i])

		ax = gca()
		if i < nspec - 1: ax.set_xticklabels([])

	subplots_adjust(hspace=0)
	xlabel("Wavelength [\AA]")

	savefig('fig_uv.eps',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()



def make_fig_uv_opt_comp(filenames):

	text_locs = [1600,3100]
	xmin = []
	xmax = []

	for i in range(len(filenames)):
		s = rd.read_spec_file(filename)

		figure(figsize=(8.3,11.6))

		nspec = len(s.spec)

		incs = ("10", "17.5", "45", "62.5", "80")

		for i in range(len(s.spec)):

			subplot(nspec,1,i+1)
			xlim(xmin[i],xmax[i])

			plot(s.wavelength, 1e12*smooth(s.spec[i] ), c="k")

			if i == 2: ylabel("Flux ($10^{-12}$~erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)")

			text(1600, 0.7*1e12*(np.max(s.spec[i])), "$i=%s^{\circ}$" % incs[i])

			ax = gca()
			if i < nspec - 1: ax.set_xticklabels([])

		subplots_adjust(hspace=0)
		xlabel("Wavelength [\AA]")

	subplots_adjust(vspace=0)
	savefig('fig_uv.eps',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()



def make_fig6():


	#files = ["sv_lk02_opt", "sv_lk02_slow_opt", "sv_lk02_slow2_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt","sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt"]
	files = ["sv_lk02_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha15_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha3_mdot3e9_rv7e10_opt.spec", "sv_lk02_alpha4_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha15_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha3_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha15_mdot3e9_rv7e11_opt.spec"]

	plt.rcParams["lines.linewidth"]= 1.5

	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", "data/atomic_macro/he_top_lines.py" ]

	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")

	nspec=9

	figure(figsize=(8.3,8.3))


	st = [r"$R_v = 7\times10^{10}$~cm", r"$R_v = 2.2\times10^{11}$~cm", r"$R_v = 7\times10^{11}$~cm"]

	for i in range(len(files)):

		s = rd.read_spec_file(files[i])

		subplot(3,3,i+1)

		print("2,3,", i+1)

		a = get_continuum(s.wavelength, s.spec[-1], lines)
		plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength), window_len=20) )

		ax = gca()

		iplot = i+1

		if iplot != 1 and iplot != 4 and iplot != 7:
			ax.set_yticklabels([])

		if iplot<7:
			ax.set_xticklabels([])

		# else:
		ax = gca()
		labels = [item.get_text() for item in ax.get_xticklabels()]

		for j in np.arange(1,len(labels),2):

			#labels[j]=""
			print(labels[j], j)

		print(labels)
		#ax.set_xticklabels(labels)

		xmin = 3410
		xmax = 4100
		xticks(np.arange(xmin , xmax, 200))

		xlim(xmin, xmax)
		ylim(0.5,1.5)
		if i<3: title(st[i])

		#if i==0:text(3500,1.2,"Balmer jump")


		# subplot(3,3,i+4)
		# a = get_continuum(s.wavelength, s.spec[-1], lines)
		# plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength)) )
		# ax = gca()
		# if i >0: ax.set_yticklabels([])
		# xlim(6513,6613)
		# ylim(0.8,3.5)

		# if i==0:text(6540,2.5,r"$H_{\alpha}$")

		# subplot(3,3,i+7)
		# a = get_continuum(s.wavelength, s.spec[-1], lines)
		# plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength)) )
		# ax = gca()
		# if i >0: ax.set_yticklabels([])
		# xlim(4616,4716)
		# ylim(0.8,3.5)

		# if i==0:text(4640,2.5,"He II $4686$")



	subplot(334)
	ylabel("Flux / Continuum", fontsize=16)
	subplot(338)
	xlabel("Wavelength [\AA]", fontsize=16)

	subplot(333)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=1.5$", fontsize=16)
	ax2.set_yticklabels([])

	subplot(336)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=3$", fontsize=16)
	ax2.set_yticklabels([])

	subplot(339)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=4$", fontsize=16)
	ax2.set_yticklabels([])

	subplots_adjust(hspace=0, wspace=0)

	savefig('fig6.eps',facecolor='w',edgecolor='w',bbox_inches='tight')
	savefig('fig6.png',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()





def vel_law (l, R_v, alpha):

	v0 = 6
	v_inf=10000

	v = (v_inf - v0) * (l / R_v)**alpha 

	v /= 1.0 + (l / R_v)**alpha 

	v += v0

	return v



def make_fig7():


	#files = ["sv_lk02_opt", "sv_lk02_slow_opt", "sv_lk02_slow2_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt","sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt"]
	files = ["sv_lk02_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha15_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha3_mdot3e9_rv7e10_opt.spec", "sv_lk02_alpha4_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha15_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha3_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha15_mdot3e9_rv7e11_opt.spec"]


	plt.rcParams["lines.linewidth"]= 1.5

	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", "data/atomic_macro/he_top_lines.py" ]

	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")

	nspec=9

	figure(figsize=(8.3,8.3))


	plt.rcParams["lines.linewidth"] = 2

	st = [r"$R_v = 7\times10^{10}$~cm", r"$R_v = 2.2\times10^{11}$~cm", r"$R_v = 7\times10^{11}$~cm"]

	for i in range(len(files)):

		s = rd.read_spec_file(files[i])

		subplot(3,3,i+1)

		print("2,3,", i+1)

		a = get_continuum(s.wavelength, s.spec[-1], lines)
		plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength), window_len=10) )

		ax = gca()

		iplot = i+1

		if iplot != 1 and iplot != 4 and iplot != 7:
			ax.set_yticklabels([])

		if iplot<7:
			ax.set_xticklabels([])

		# else:
		ax = gca()
		labels = [item.get_text() for item in ax.get_xticklabels()]

		for j in np.arange(1,len(labels),2):

			#labels[j]=""
			print (labels[j], j)

		print (labels)
		#ax.set_xticklabels(labels)

		xmin = 6513
		xmax = 6613
		xticks(np.arange(xmin , xmax, 20))

		xlim(xmin, xmax)
		ylim(0.5,10)
		if i<3: title(st[i])



	subplot(334)
	ylabel("Flux / Continuum", fontsize=16)
	subplot(338)
	xlabel("Wavelength [\AA]", fontsize=16)

	subplot(333)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=1.5$", fontsize=16)
	ax2.set_yticklabels([])

	subplot(336)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=3$", fontsize=16)
	ax2.set_yticklabels([])

	subplot(339)
	ax2 = gca().twinx()
	ax2.set_ylabel(r"$\alpha=4$", fontsize=16)
	ax2.set_yticklabels([])

	subplots_adjust(hspace=0, wspace=0)

	savefig('fig7.eps',facecolor='w',edgecolor='w',bbox_inches='tight')
	savefig('fig7.png',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()

def demo_images_side_by_side(ax):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)

    #Z, extent = get_demo_image()
    ax2 = divider.new_horizontal(size="100%", pad=0.05)
    fig1 = ax.get_figure()
    fig1.add_axes(ax2)

    return ax, ax2
    # ax.imshow(Z, extent=extent, interpolation="nearest")
    # ax2.imshow(Z, extent=extent, interpolation="nearest")
    # for tl in ax2.get_yticklabels():
    #     tl.set_visible(False)


def make_fig8(files, nsubfigs, nx = 3, ny = 3, xmin = 6513, xmax = 6613, size = (16,8), savename="fig8", ymin=0.5, ymax=7.9, dtick=50):

	#files = ["sv_lk02_opt", "sv_lk02_slow_opt", "sv_lk02_slow2_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt","sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt"]
	#_files = ["sv_lk02_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha15_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha3_mdot3e9_rv7e10_opt.spec", "sv_lk02_alpha4_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha15_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha3_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha15_mdot3e9_rv7e11_opt.spec"]

	plt.rcParams["lines.linewidth"]= 1.5
	plt.rcParams["text.usetex"]="True"

	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", 
	                     "data/atomic_macro/he_top_lines.py" ]
	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")


	nspec = nx * ny

	if len(files) != nsubfigs * nspec:
		print("Error, bad dimensions!")
		return (-1)


	figure(figsize=size)
	


	# some figure parameters controlling spacing
	fig_gap = 0.03
	margin = 0.05
	fig_size = (1.0 - fig_gap*(nsubfigs - 1) - margin*2) / nsubfigs


	glist = []
	# create the gridspec instances
	for i in range(nsubfigs):

		l = margin + (fig_gap * i) + (fig_size*i)
		r = margin + (fig_gap * i) + (fig_size*(i+1))

		glist.append(gridspec.GridSpec(nx, ny))
		glist[i].update(left = l, right = r, top = 0.85, wspace=0.0, hspace = 0.0)


	mdots = [r'$\dot{m}_W=10^{-9}$', r'$\dot{m}_W=3\times10^{-9}$', r'$\dot{m}_W=3\times10^{-9}$']
	st = [r"$7\times10^{10}$~cm", r"$2.2\times10^{11}$~cm", r"$7\times10^{11}$~cm"]

	pos_h = ["l","m","r","l","m","r","l","m","r"]
	pos_v = ["t","t","t","m","m","m","b","b","b"]


	for n in range(nsubfigs):

		istart = n * nspec 

		for i in range(nspec):

			s = rd.read_spec_file(files[istart + i])
			subplot(glist[n][i])


			print (i)

			if pos_v[i] == "b" and pos_h[i]=="m":
				xlabel("Wavelength [\AA]", fontsize=20)

			if n==0 and pos_h[i]=="l" and pos_v[i]=="m":
				ylabel("Flux / Continuum", fontsize=20)

			a = get_continuum(s.wavelength, s.spec[-1], lines)

			noplot = (n == 2 and i == 8)

			if noplot == False:
				plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength), 
					 window_len=10) )


			ax = gca()
			iplot = i+1

			if pos_h[i]!="l":
				ax.set_yticklabels([])

			if pos_v[i]!="b":
				ax.set_xticklabels([])

			if pos_h[i]=="m" and pos_v[i]=="t":
				ax = gca()
				text(0.5, 1.3,mdots[n], horizontalalignment='center',
					 verticalalignment='center',transform=ax.transAxes, 
					 fontsize = 16)

			# else:
			labels = [item.get_text() for item in ax.get_xticklabels()]

			#ax.set_xticklabels(labels)

			#xmin = 6513
			#xmax = 6613

			xticks(np.arange(xmin , xmax, dtick))

			xlim(xmin, xmax)
			ylim(ymin,ymax)
			#yticks(np.arange(0.5,7.9,1))
			if i<3: title(st[i])

			if n == 2 and pos_h[i] =="r":
				ax2 = gca().twinx()

				if i == 2:
					ax2.set_ylabel(r"$\alpha=1.5$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 5:
					ax2.set_ylabel(r"$\alpha=3$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 8:
					ax2.set_ylabel(r"$\alpha=4$", fontsize=16)
					ax2.set_yticklabels([])


	savefig(savename+'.eps',facecolor='w',edgecolor='w',bbox_inches='tight')
	savefig(savename+'.png')


	show()

	return 0





def make_fig9(files, nsubfigs, nx = 3, ny = 3, xmin = 6513, xmax = 6613, size = (16,8)):

	#files = ["sv_lk02_opt", "sv_lk02_slow_opt", "sv_lk02_slow2_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt","sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt"]
	#_files = ["sv_lk02_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha15_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha3_mdot3e9_rv7e10_opt.spec", "sv_lk02_alpha4_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha15_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha3_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha15_mdot3e9_rv7e11_opt.spec"]

	plt.rcParams["lines.linewidth"]= 1.5
	plt.rcParams["text.usetex"]="True"

	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", 
	                     "data/atomic_macro/he_top_lines.py" ]
	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")


	nspec = nx * ny

	if len(files) != nsubfigs * nspec:
		print("Error, bad dimensions!")
		return (-1)


	figure(figsize=size) 
	


	# some figure parameters controlling spacing
	fig_gap = 0.03
	margin = 0.05
	fig_size = (1.0 - fig_gap*(nsubfigs - 1) - margin*2) / nsubfigs


	glist = []
	# create the gridspec instances
	for i in range(nsubfigs):

		l = margin + (fig_gap * i) + (fig_size*i)
		r = margin + (fig_gap * i) + (fig_size*(i+1))

		glist.append(gridspec.GridSpec(nx, ny))
		glist[i].update(left = l, right = r, top = 0.85, wspace=0.0, hspace = 0.0)


	mdots = [r'$\dot{m}_W=10^{-9}$', r'$\dot{m}_W=3\times10^{-9}$', r'$\dot{m}_W=3\times10^{-9}$']
	st = [r"$7\times10^{10}$~cm", r"$2.2\times10^{11}$~cm", r"$7\times10^{11}$~cm"]

	pos_h = ["l","m","r","l","m","r","l","m","r"]
	pos_v = ["t","t","t","m","m","m","b","b","b"]


	for n in range(nsubfigs):

		istart = n * nspec 

		for i in range(nspec):

			s = rd.read_spec_file(files[istart + i])
			subplot(glist[n][i])


			print(i)

			if pos_v[i] == "b" and pos_h[i]=="m":
				xlabel("Wavelength [\AA]", fontsize=20)

			if n==0 and pos_h[i]=="l" and pos_v[i]=="m":
				ylabel("Flux / Continuum", fontsize=20)

			a = get_continuum(s.wavelength, s.spec[-1], lines)

			noplot = (n == 2 and i == 8)

			if noplot == False:
				plot(s.wavelength, smooth(s.spec[-1]/make_f(a, s.wavelength), 
					 window_len=10) )


			ax = gca()
			iplot = i+1

			if pos_h[i]!="l":
				ax.set_yticklabels([])

			if pos_v[i]!="b":
				ax.set_xticklabels([])

			if pos_h[i]=="m" and pos_v[i]=="t":
				ax = gca()
				text(0.5, 1.3,mdots[n], horizontalalignment='center',
					 verticalalignment='center',transform=ax.transAxes, 
					 fontsize = 16)

			# else:
			labels = [item.get_text() for item in ax.get_xticklabels()]

			#ax.set_xticklabels(labels)

			xmin = 3400
			xmax = 4000
			xticks(np.arange(xmin , xmax, 200))

			xlim(xmin, xmax)
			ylim(0.2,2)
			#yticks(np.arange(0.5,7.9,1))
			if i<3: title(st[i])

			if n == 1 and pos_h[i] =="r":
				ax2 = gca().twinx()

				if i == 2:
					ax2.set_ylabel(r"$\alpha=1.5$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 5:
					ax2.set_ylabel(r"$\alpha=3$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 8:
					ax2.set_ylabel(r"$\alpha=4$", fontsize=16)
					ax2.set_yticklabels([])


	savefig('fig9.eps',facecolor='w',edgecolor='w',bbox_inches='tight')
	savefig('fig9.png',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()

	return 0







def make_fig_uv(files, nsubfigs, nx = 3, ny = 3, xmin = 1500, xmax = 1600, size = (16,8)):

	#files = ["sv_lk02_opt", "sv_lk02_slow_opt", "sv_lk02_slow2_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt","sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha3_rtimes3_mdot_opt"]
	#_files = ["sv_lk02_opt", "sv_lk02_alpha3_opt", "sv_lk02_alpha3_rtimes3_opt", "sv_lk02_alpha15_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha3_mdot3e9_rv7e10_opt.spec", "sv_lk02_alpha4_mdot3e9_rv7e10_opt.spec","sv_lk02_alpha15_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha3_mdot3e9_rv2e11_opt.spec","sv_lk02_alpha15_mdot3e9_rv7e11_opt.spec"]

	plt.rcParams["lines.linewidth"]= 1.5
	plt.rcParams["text.usetex"]="True"

	atomic_filenames = [ "data/atomic_macro/h10_lines.py", "data/atomic73/lines_linked_ver_2.py", 
	                     "data/atomic_macro/he_top_lines.py" ]
	lines = get_linelist(atomic_filenames)

	incs = ("10", "17.5", "45", "62.5", "80")


	nspec = nx * ny


	figure(figsize=size)
	


	# some figure parameters controlling spacing
	fig_gap = 0.03
	margin = 0.05
	fig_size = (1.0 - fig_gap*(nsubfigs - 1) - margin*2) / nsubfigs


	glist = []
	# create the gridspec instances
	for i in range(nsubfigs):

		l = margin + (fig_gap * i) + (fig_size*i)
		r = margin + (fig_gap * i) + (fig_size*(i+1))

		glist.append(gridspec.GridSpec(nx, ny))
		glist[i].update(left = l, right = r, top = 0.85, wspace=0.0, hspace = 0.0)


	mdots = [r'$\dot{m}_W=2\times10^{-9}, i=45$', r'$\dot{m}_W=2\times10^{-9}, i=62.5$', r'$\dot{m}_W=2\times10^{-9}, i=80$']
	st = [r"$7\times10^{10}$~cm", r"$2.2\times10^{11}$~cm", r"$7\times10^{11}$~cm"]


	pos_h = ["l","m","r","l","m","r","l","m","r"]
	pos_v = ["t","t","t","m","m","m","b","b","b"]

	angle = [-3,-2, -1]

	for n in range(nsubfigs):

		istart = 0



		for i in range(nspec):

			s = rd.read_spec_file(files[istart + i])
			subplot(glist[n][i])


			print(i)

			if pos_v[i] == "b" and pos_h[i]=="m":
				xlabel("Wavelength [\AA]", fontsize=20)

			if n==0 and pos_h[i]=="l" and pos_v[i]=="m":
				ylabel("Flux / Continuum", fontsize=20)

			#a = get_continuum(s.wavelength, s.spec[-1], lines, lmin=1300,lmax=1800)

			#noplot = (n == 2 and i == 8)

			#if noplot == False:
			plot(s.wavelength, smooth(s.spec[angle[n]]/s.spec[angle[n]][4000], window_len=10) )

			
			ax = gca()
			iplot = i+1

			if pos_h[i]!="l":
				ax.set_yticklabels([])

			if pos_v[i]!="b":
				ax.set_xticklabels([])

			if pos_h[i]=="m" and pos_v[i]=="t":
				ax = gca()
				text(0.5, 1.3,mdots[n], horizontalalignment='center',
					 verticalalignment='center',transform=ax.transAxes, 
					 fontsize = 16)

			# else:
			labels = [item.get_text() for item in ax.get_xticklabels()]

			#ax.set_xticklabels(labels)

			#xmin = 3400
			#xmax = 4000
			xticks(np.arange(xmin , xmax, 50))

			xlim(xmin, xmax)
			ylim(0,10)
			#yticks(np.arange(0.5,7.9,1))
			if i<3: title(st[i])

			if n == 2 and pos_h[i] =="r":
				ax2 = gca().twinx()

				if i == 2:
					ax2.set_ylabel(r"$\alpha=1.5$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 5:
					ax2.set_ylabel(r"$\alpha=3$", fontsize=16)
					ax2.set_yticklabels([])

				if i == 8:
					ax2.set_ylabel(r"$\alpha=4$", fontsize=16)
					ax2.set_yticklabels([])
			

	savefig('figuv.eps',facecolor='w',edgecolor='w',bbox_inches='tight')
	savefig('figuv.png',facecolor='w',edgecolor='w',bbox_inches='tight')

	show()

	return s








