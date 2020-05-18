"""
Copyright (C) 2020 Gebri Mishtaku

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses.
"""

import os
import json
import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve
from utils import *



class CenterFinder(object):

	def __init__(self, galaxy_file: str, kernel_radius: float, vote_threshold: float, 
		g_wtd: bool = False, params_file: str = None, save: bool = False, printout: bool = False):

		self.kernel_radius = kernel_radius
		self.vote_threshold = vote_threshold

		self.filename = '.'.join(galaxy_file.split('.')[:-1])
		self.save = save
		if self.save:
			self.savename = f'out_{self.filename}/'
		self.printout = printout

		# loads galaxy data arrays
		if not g_wtd:
			self.G_ra, self.G_dec, self.G_redshift = load_data(galaxy_file)
			self.G_weights = np.ones(len(self.G_ra), dtype=float)
		else:
			self.G_ra, self.G_dec, self.G_redshift, self.G_weights = load_data_weighted(galaxy_file)

		# gets cosmology and other hyperparameters
		self.cosmology, self.grid_spacing = load_hyperparameters(params_file)
		# calculates lookup tables for fast conversion from r to z and vice versa
		self.LUT_radii, self.LUT_redshifts = interpolate_r_z(self.G_redshift.min(), self.G_redshift.max(), self.cosmology)

		self.G_radii = self.LUT_radii(self.G_redshift)



	def _kernel(self, additional_thickness: float = 0., show_kernel: bool = False) -> np.ndarray:
		# this is the number of bins in each dimension axis
		# this calculation ensures an odd numbered gridding
		# the kernel construction has a distinct central bin on any given run
		kernel_bin_count = int(2 * np.ceil(self.kernel_radius / self.grid_spacing) + 1)

		# this is the kernel inscribed radius in index units
		inscribed_r_idx_units = self.kernel_radius / self.grid_spacing
		inscribed_r_idx_units_upper_bound = inscribed_r_idx_units + 0.5 + additional_thickness
		inscribed_r_idx_units_lower_bound = inscribed_r_idx_units - 0.5 - additional_thickness

		# central bin index, since the kernel is a cube this can just be one int
		kernel_center_index = int(kernel_bin_count / 2)
		kernel_center = np.array([kernel_center_index, ] * 3)

		# this is where the magic happens: each bin at a radial distance of bao_radius from the
		# kernel's center gets assigned a 1 and all other bins get a 0
		kernel_grid = np.array([[[1 if (np.linalg.norm(np.array([i, j, k]) - kernel_center) >= inscribed_r_idx_units_lower_bound
										and np.linalg.norm(np.array([i, j, k]) - kernel_center) < inscribed_r_idx_units_upper_bound)
								  else 0
								  for k in range(kernel_bin_count)]
								 for j in range(kernel_bin_count)]
								for i in range(kernel_bin_count)])

		if self.printout:
			print('Kernel constructed successfully...')
			print('Number of kernel bins that contain surface:', len(kernel_grid[kernel_grid == 1]))
			print('Number of empty kernel bins:', len(kernel_grid[kernel_grid == 0]))

		# this is here for future sanity checks, it shows the kernel in 3d
		# with blue disks in kernel bins containing spherical surface
		if show_kernel:
			color = 'cornflowerblue'
			fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
			ax.scatter(*np.where(kernel_grid == 1), c=color)
			plt.show()

		return kernel_grid



	def _vote(self, background_subtract: bool = False, plot: bool = False):

		xyzs = sky2cartesian(self.G_ra, self.G_dec, self.G_redshift, self.LUT_radii) # galaxy x, y and z coordinates
		self.galaxies_cartesian = np.array(xyzs).T  # each galaxy is represented by (x, y, z)

		# gets the 3d histogram (density_grid) and the grid bin coordintes in cartesian (grid_edges)
		bin_counts_3d = np.array([np.ceil((xyzs[i].max() - xyzs[i].min()) / self.grid_spacing) 
									for i in range(len(xyzs))], dtype=int)
		density_grid, self.density_grid_edges = np.histogramdd(self.galaxies_cartesian, bins=bin_counts_3d)

		if self.printout:
			print('Histogramming completed successfully...')
			print('Density grid shape:', density_grid.shape)

		# subtracts the background
		if background_subtract:
			background, _ = project_and_sample(density_grid, self.density_grid_edges, printout=printout)
			density_grid -= background
			density_grid[density_grid < 0.] = 0.

		# makes the kernel for scanning over the density grid
		kernel_grid = self._kernel()

		# this scans the kernel over the whole volume of the galaxy density grid
		# calculates the tensor inner product of the two at each step
		# and finally stores this value as the number of voters per that bin in the observed grid
		self.centers_grid = np.round(fftconvolve(density_grid, kernel_grid, mode='same'))
		
		if self.printout:
			print('Observed grid shape:', self.centers_grid.shape)
			print('Maximum number of voters per single bin:', self.centers_grid.max())
			print('Minimum number of voters per single bin:', self.centers_grid.min())

		# save whole grid without a vote cut
		if self.save:
			np.save(self.savename + f'centers_grid_r_{self.kernel_radius}_no_cut.npy', 
				self.centers_grid)



	def find_centers(self):

		self._vote()

		# applies vote threshold and saves the found centers list as .fits catalog
		centers_indices = np.asarray(self.centers_grid >= self.vote_threshold).nonzero()
		self.C_weights = self.centers_grid[centers_indices]
		delattr(self, 'centers_grid')
		centers_bin_coords = np.array([(self.density_grid_edges[i][:-1] + self.density_grid_edges[i][1:]) / 2 for i in range(len(self.density_grid_edges))])
		delattr(self, 'density_grid_edges')
		C_xyzs = np.array([centers_bin_coords[i][centers_indices[i]] for i in range(len(centers_indices))])
		self.C_ra, self.C_dec, self.C_redshift, _ = cartesian2sky(*C_xyzs, self.LUT_redshifts)
		
		if self.save:
			savename = self.savename + f'found_centers_r_{self.kernel_radius}_cut_{self.vote_threshold}.fits'
			save_data_weighted(savename, self.C_ra, self.C_dec, self.C_redshift, self.C_weights)







