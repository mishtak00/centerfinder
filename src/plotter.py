import numpy as np
import matplotlib.pyplot as plt
plt.rc('figure', figsize=[9,9])
plt.rc('font', family='serif')
plt.rc('axes', titlesize=18)
plt.rc('axes', labelsize=12)
plt.rc('xtick', top=True)
plt.rc('xtick.minor', visible=True)
plt.rc('ytick', right=True)
plt.rc('ytick.minor', visible=True)



def _plot_slice(cf, bounds):

	# user specified bounds on all coords
	if len(bounds)==6:
		zlow, zhi, ralow, rahi, declow, dechi = bounds
	# specified bounds only on redshift, max bounds on ra and dec
	elif len(bounds)==2:
		zlow, zhi = bounds
		ralow, rahi = cf.G_ra.min(), cf.G_ra.max()
		declow, dechi = cf.G_dec.min(), cf.G_dec.max()
	# specified no bounds, max bounds assumed on all coords
	elif len(bounds)==0:
		zlow, zhi = cf.G_redshift.min(), cf.G_redshift.max()
		ralow, rahi = cf.G_ra.min(), cf.G_ra.max()
		declow, dechi = cf.G_dec.min(), cf.G_dec.max()
	else:
		raise ValueError('Plotting Error: Number of plotting function arguments '\
			f'= {len(bounds)}, but expected either 0, 2 or 6')

	if zlow>zhi or ralow>rahi or declow>dechi:
		raise ValueError('Plotting Error: Lower bound of a coordinate '\
			'cannot be greater than higher bound of said coordinate.')

	# cannot plot 0 centers
	if len(cf.C_weights)==0:
		raise ValueError('Plotting Error: Centers list is empty.')

	fig, ax = plt.subplots()

	# gets centers to plot
	slice_indices = np.asarray((cf.C_redshift>=zlow) & (cf.C_redshift<=zhi)).nonzero()
	ra_slice, dec_slice = cf.C_ra[slice_indices], cf.C_dec[slice_indices]
	wts_slice = cf.C_weights[slice_indices]
	h, ra, dec, _ = ax.hist2d(ra_slice, dec_slice, weights=wts_slice, bins=100, 
		range=[[ralow,rahi],[declow,dechi]], alpha=0.7)

	# gets galaxies to plot
	slice_indices = np.asarray((cf.G_redshift>=zlow) & (cf.G_redshift<=zhi) 
		& (cf.G_ra>=ralow) & (cf.G_ra<=rahi)
		& (cf.G_dec>=declow) & (cf.G_dec<=dechi)).nonzero()
	ra_slice, dec_slice = cf.G_ra[slice_indices], cf.G_dec[slice_indices]
	ax.plot(ra_slice, dec_slice, 'r*', ms=4)

	# plots circles around hihgly voted centers
	cell_density_floor = .8*h.max()
	circ_indices = np.asarray(h>=cell_density_floor).nonzero()
	# bounds the number of circles on screen
	while len(circ_indices[0])>20:
		cell_density_floor*=1.01
		circ_indices = np.asarray(h>=cell_density_floor).nonzero()
	ra, dec = ra[circ_indices[0]], dec[circ_indices[1]]
	circ_centers = list(zip(ra, dec))
	for xy in circ_centers:
		# TODO: automatize calculation for angle subtended by R0 given z
		ax.add_artist(plt.Circle(xy, 4.5, ls='--', fill=False, alpha=.5))

	plt.title('Plot of centers and galaxies between $z={:.2f}$ and '\
		'$z={:.2f}$'.format(zlow, zhi))
	plt.xlabel('$RA$ $[\\degree]$')
	plt.ylabel('$DEC$ $[\\degree]$')

	if cf.save:
		savename = cf.savename + 'slice_plot_z_{:.2f}_{:.2f}.png'.format(zlow,zhi)
		plt.savefig(savename, dpi=300)

	plt.show()


def _plot_coord_hist(cf, which: str):

	if which=='RA':
		gcoord, ccoord = cf.G_ra, cf.C_ra
		plt.xlabel(which + r' $[\degree]$')
	elif which=='DEC':
		gcoord, ccoord = cf.G_dec, cf.C_dec
		plt.xlabel(which + r' $[\degree]$')
	elif which=='Z':
		gcoord, ccoord = cf.G_redshift, cf.C_redshift
		plt.xlabel(which + r' $[h^{-1}Mpc]$')
	elif which=='R':
		gcoord, ccoord = cf.G_radii, cf.C_radii
		plt.xlabel(which + r' $[h^{-1}Mpc]$')
	else:
		raise ValueError('Plotting Error: Can only plot histogram for one of '\
		 'RA, DEC, Z, R; instead found "{}"'.format(which))

	plt.hist(gcoord, label='Galaxies ' + which)
	plt.hist(ccoord, label='Centers ' + which, histtype='step')
	plt.title('Galaxy and center unweighted counts by ' + which)
	plt.ylabel('Count')
	plt.legend()

	if cf.save:
		savename = cf.savename + '{}_hist_r_{}_cut_{}.png'.format(which,
			cf.kernel_radius,cf.vote_threshold)
		plt.savefig(savename, dpi=300)

	plt.show()




