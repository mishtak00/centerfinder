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
from math import inf
from astropy.io import fits
from scipy.signal import fftconvolve
from .utils import *
from .kernel import Kernel



class CenterFinder():


	def __init__(self, galaxy_file: str, wtd: bool, params_file: str, save: bool, printout: bool,
		kernel_radius: float = 110., kernel_type: str = 'step', 
		kernel_args: list = [], vote_threshold: float = -inf):

		self.kernel_radius = kernel_radius
		self.kernel_type = kernel_type
		self.kernel_args = kernel_args
		self.show_kernel = False
		self.vote_threshold = vote_threshold

		self.filename = '.'.join(galaxy_file.split('.')[:-1])
		self.save = save
		self.savename = f'out_{self.filename}/'
		self.printout = printout

		# loads galaxy data arrays
		if not wtd:
			self.G_ra, self.G_dec, self.G_redshift = load_data(galaxy_file)
			self.G_weights = np.ones(len(self.G_ra), dtype=float)
		else:
			self.G_ra, self.G_dec, self.G_redshift, self.G_weights = load_data_weighted(galaxy_file)

		# gets cosmology and other hyperparameters
		self.cosmology, self.grid_spacing = load_hyperparameters(params_file)
		# calculates lookup tables for fast conversion from r to z and vice versa
		self.LUT_radii, self.LUT_redshifts = interpolate_r_z(self.G_redshift.min(), 
			self.G_redshift.max(), self.cosmology)

		self.G_radii = self.LUT_radii(self.G_redshift)



	def __str__(self):
		return 'CenterFinder object\n'\
				f'Galaxy data file: {self.filename}\n'\
				f'Kernel radius: {self.kernel_radius}\n'\
				f'Vote threshold: {self.vote_threshold}\n'\
				f'RA range: [{self.G_ra.min()}, {self.G_ra.max()}]\n'\
				f'DEC range: [{self.G_dec.min()}, {self.G_dec.max()}]\n'\
				f'Z range: [{self.G_redshift.min()}, {self.G_redshift.max()}]'



	def set_kernel_radius(self, kr: float):
		self.kernel_radius = kr


	def set_kernel_type(self, kt: str, args):
		self.kernel_type = kt
		self.kernel_args = args


	def set_show_kernel(self, sk: bool):
		self.show_kernel = sk


	def set_vote_threshold(self, vt: float):
		self.vote_threshold = vt



	def _vote(self, dencon: tuple = (False, False), overden: bool = False, 
		premade_kernel: np.ndarray = None):

		xyzs = sky2cartesian(self.G_ra, self.G_dec, self.G_redshift, self.LUT_radii) # galaxy x, y and z coordinates
		self.galaxies_cartesian = np.array(xyzs).T  # each galaxy is represented by (x, y, z)

		# gets the 3d histogram (density_grid) and the grid bin coordintes in cartesian (grid_edges)
		bin_counts_3d = np.array([np.ceil((xyzs[i].max() - xyzs[i].min()) / self.grid_spacing) 
									for i in range(len(xyzs))], dtype=int)

		# histograms the data points in real space with given weights
		density_grid, self.density_grid_edges = np.histogramdd(self.galaxies_cartesian, 
			bins=bin_counts_3d, weights=self.G_weights)

		if self.printout:
			print('Histogramming completed successfully...')
			print('Density grid shape:', density_grid.shape)

		# makes expected grid and subtracts it from the density grid
		if dencon[0]:
			background = self._project_and_sample(density_grid, self.density_grid_edges)
			density_grid -= background
			del background
			# keep or discard negative valued weights
			# dencon[1] set to True means keep negative weights
			if not dencon[1]:
				if self.printout:
					print('Discarding all negative weights in density grid...')
				density_grid[density_grid < 0.] = 0.
			if self.printout:
				print('Background subtraction completed successfully...')
		
		# calculates avg density of all nonempty grid cells of the 
		# weighted density field and subtracts it from the density field
		elif overden:
			denavg = np.average(density_grid[density_grid!=0])
			density_grid[density_grid!=0] -= denavg
			density_grid[density_grid!=0] /= denavg
			if self.printout:
				print('Overdensity calculation completed successfully...')

		if self.printout:
			print('Minimum and maximum values of density field grid cells: '\
				'[{}, {}]'.format(density_grid.min(), density_grid.max()))

		# sets kernel for the refinement process
		if premade_kernel is None:
			# makes the kernel for scanning over the density grid
			kernel = Kernel(self.kernel_type, self.kernel_radius, self.grid_spacing,
				self.printout, self.show_kernel, *self.kernel_args)
			# for refinement process, the following allows kernel jittering
			# box around each center is 5% larger than the kernel outer radius
			self.centerbox_r_upper_bound_idx_units = np.round(kernel.\
				kernel_r_idx_units_upper_bound*1.05)
		else:
			kernel = premade_kernel

		# TODO: for refinement DEV ONLY
		self.density_grid = density_grid

		# this scans the kernel over the whole volume of the galaxy density grid
		# calculates the tensor inner product of the two at each step
		# and finally stores this value as the number of voters per that bin in the centers grid
		self.centers_grid = fftconvolve(density_grid, kernel.get_grid(), mode='same')
		
		if self.printout:
			print('Voting procedure completed successfully...')
			print('Centers grid shape:', self.centers_grid.shape)
			print('Maximum number of votes per single bin:', self.centers_grid.max())
			print('Minimum number of votes per single bin:', self.centers_grid.min())

		# save whole grid without a vote cut
		if self.save:
			np.save(self.savename + f'centers_grid_r_{self.kernel_radius}_no_cut.npy', 
				self.centers_grid)



	def _project_and_sample(self, grid: np.ndarray, grid_edges: list) -> (np.ndarray, tuple):
		
		# TODO: these are unnecessary, remove and reference self's attribute
		bin_centers_edges_xs, bin_centers_edges_ys, bin_centers_edges_zs = \
			np.array([(grid_edges[i][:-1] + grid_edges[i][1:]) / 2 for i in range(len(grid_edges))])

		# TODO: remove, unnecessary
		# if self.save:
		# 	np.save(self.savename + '_xbins.npy', bin_centers_edges_xs)
		# 	np.save(self.savename + '_ybins.npy', bin_centers_edges_ys)
		# 	np.save(self.savename + '_zbins.npy', bin_centers_edges_zs)

		bin_centers_xs, bin_centers_ys, bin_centers_zs = np.array([(x, y, z) 
			for x in bin_centers_edges_xs 
			for y in bin_centers_edges_ys 
			for z in bin_centers_edges_zs]).T
		del bin_centers_edges_xs, bin_centers_edges_ys, bin_centers_edges_zs
		if self.printout:
			print('Number of bin centers in cartesian coordinates:', len(bin_centers_xs))
		"""
		Why can we be sure that it is okay to interpolate the radii 
		and redshift values for these bin centers coordinates?
		Because we know that the range of values of the bin centers 
		is exactly in between the min and the max of the grid bin 
		edges x, y, z. The radii come from the 3d euclidian distance, 
		which preserves this relationship (convex function of x,y,z), 
		and thus it is fine to use the beforehand-calculated interpolation 
		lookup table to find the redshifts from the radii.
		"""
		bin_centers_ra, bin_centers_dec, _, bin_centers_radii = \
			cartesian2sky(bin_centers_xs, 
							bin_centers_ys, 
							bin_centers_zs, 
							self.LUT_redshifts, 
							self.G_ra.min(), 
							self.G_ra.max())
		del bin_centers_xs, bin_centers_ys, bin_centers_zs
		if self.printout:
			print('Number of bin centers in sky coordinates:', len(bin_centers_ra))

		# total number of votes
		N_tot = np.sum(grid)
		if self.printout:
			print('Total number of votes:', N_tot)

		# get volume adjustment grid, the differentials in sky coordinate dimensions and the number of bins in each dimension
		vol_adjust_ratio_grid, d_r, d_alpha, d_delta, N_bins_r, N_bins_alpha, N_bins_delta = \
			self._volume_adjustment(bin_centers_radii, bin_centers_ra, bin_centers_dec, grid.shape)

		# alpha-delta and z counts
		N_bins_x, N_bins_y, N_bins_z = grid.shape[0], grid.shape[1], grid.shape[2]
		sky_coords_grid_shape = (N_bins_x, N_bins_y, N_bins_z, 3)  # need to store a triple at each grid bin
		sky_coords_grid = np.array(list(zip(bin_centers_ra, bin_centers_dec, 
			bin_centers_radii))).reshape(sky_coords_grid_shape)
		if self.printout:
			print('Shape of grid containing sky coordinates of observed grid\'s bin centers:', 
				sky_coords_grid.shape)

		# getting some variables ready for the projection step
		alpha_min = bin_centers_ra.min()
		d_alpha = np.rad2deg(d_alpha)
		delta_min = bin_centers_dec.min()
		d_delta = np.rad2deg(d_delta)
		r_min = bin_centers_radii.min()
		del bin_centers_ra, bin_centers_dec, bin_centers_radii

		# vectorial computation of the sky indices
		sky_coords_grid[:, :, :, 0] = (sky_coords_grid[:, :, :, 0] - alpha_min) // d_alpha
		sky_coords_grid[:, :, :, 1] = (sky_coords_grid[:, :, :, 1] - delta_min) // d_delta
		sky_coords_grid[:, :, :, 2] = (sky_coords_grid[:, :, :, 2] - r_min) // d_r
		sky_coords_grid = sky_coords_grid.astype(int)

		# TODO: the condition here should be >= rather than ==
		# the following fixes any indices that lie beyond the outer 
		# walls of the sky grid by pulling them to the wall
		sky_coords_grid[:, :, :, 0][sky_coords_grid[:, :, :, 0] == N_bins_alpha] = N_bins_alpha - 1
		sky_coords_grid[:, :, :, 1][sky_coords_grid[:, :, :, 1] == N_bins_delta] = N_bins_delta - 1
		sky_coords_grid[:, :, :, 2][sky_coords_grid[:, :, :, 2] == N_bins_r] = N_bins_r - 1

		alpha_delta_grid, r_grid = self._alpha_delta_r_projections_from_grid(grid, 
			N_bins_x, N_bins_y, N_bins_z, sky_coords_grid, N_bins_alpha, N_bins_delta, N_bins_r)
		if self.printout:
			print('Shape of alpha-delta grid:', alpha_delta_grid.shape)
			print('Shape of r grid:', r_grid.shape)
			print('Maximum number of voters per single bin in alpha-delta grid:', alpha_delta_grid.max())
			print('Minimum number of voters per single bin in alpha-delta grid:', alpha_delta_grid.min())
			print('Maximum number of voters per single bin in r grid:', r_grid.max())
			print('Minimum number of voters per single bin in r grid:', r_grid.min())
			print('N_tot_observed = N_tot_alpha_delta = N_tot_r:', 
				N_tot == np.sum(alpha_delta_grid) == np.sum(r_grid))

		# TODO: meshgrid this
		expected_grid = np.array([[[alpha_delta_grid[sky_coords_grid[i, j, k, 0], sky_coords_grid[i, j, k, 1]]
									 * r_grid[sky_coords_grid[i, j, k, 2]]
									for k in range(N_bins_z)]
								   for j in range(N_bins_y)]
								  for i in range(N_bins_x)])

		expected_grid /= N_tot  # normalization
		expected_grid *= vol_adjust_ratio_grid  # volume adjustment
		if self.printout:
			print('Expected grid shape:', expected_grid.shape)
			print('Maximum number of expected votes:', expected_grid.max())
			print('Minimum number of expected votes:', expected_grid.min())

		# if self.save:
		# 	np.save(self.savename + "_exp_grid.npy", expected_grid)

		return expected_grid



	def _volume_adjustment(self, bin_centers_radii: np.array, bin_centers_ra: np.array, 
		bin_centers_dec: np.array, observed_grid_shape: tuple) -> np.ndarray:
		
		# radius
		mid_r = (bin_centers_radii.max() + bin_centers_radii.min()) / 2
		delta_r = bin_centers_radii.max() - bin_centers_radii.min()
		N_bins_r = int(np.ceil(delta_r / self.grid_spacing))
		d_r = self.grid_spacing
		r_sqr = bin_centers_radii ** 2

		# alpha
		delta_alpha = np.deg2rad(bin_centers_ra.max() - bin_centers_ra.min())
		N_bins_alpha = int(np.ceil((delta_alpha * mid_r) / self.grid_spacing))
		d_alpha = delta_alpha / N_bins_alpha

		# delta
		delta_delta = np.deg2rad(bin_centers_dec.max() - bin_centers_dec.min())
		N_bins_delta = int(np.ceil((delta_delta * mid_r) / self.grid_spacing))
		d_delta = delta_delta / N_bins_delta
		cos_delta = np.cos(np.deg2rad(bin_centers_dec))

		# angular volume differential
		dV_ang = d_alpha * cos_delta * d_delta * r_sqr * d_r
		# euclidean volume differential
		dV_xyz = self.grid_spacing ** 3
		# volume adjustment ratio grid; contains the volume adjustment ratio per each bin in the expected grid
		vol_adjust_ratio_grid = (dV_xyz / dV_ang).reshape(observed_grid_shape)

		if self.printout:
			print('Number of bins in r:', N_bins_r)
			print('Number of bins in alpha:', N_bins_alpha)
			print('Number of bins in delta:', N_bins_delta)
			print('Volume adjustment ratio grid shape:', vol_adjust_ratio_grid.shape)

		return vol_adjust_ratio_grid, d_r, d_alpha, d_delta, N_bins_r, N_bins_alpha, N_bins_delta



	def _alpha_delta_r_projections_from_grid(self, grid: np.ndarray, N_bins_x: int, N_bins_y: int, N_bins_z: int, 
		sky_coords_grid: np.ndarray, N_bins_alpha: int, N_bins_delta: int, N_bins_r: int) -> (np.ndarray, np.ndarray):

		alpha_delta_grid = np.zeros((N_bins_alpha, N_bins_delta))
		r_grid = np.zeros((N_bins_r,))
		for i in range(N_bins_x):
			for j in range(N_bins_y):
				for k in range(N_bins_z):
					alpha_delta_grid[sky_coords_grid[i, j, k, 0], sky_coords_grid[i, j, k, 1]] += grid[i, j, k]
					r_grid[sky_coords_grid[i, j, k, 2]] += grid[i, j, k]

		return alpha_delta_grid, r_grid



	def find_centers(self, dencon: bool, overden: bool):
		"""
		Identifies BAO centers by applying the voting procedure.
		Applies vote threshold and saves the found centers list as .fits catalog.
		"""

		if self.printout:
			print(self)

		self._vote(dencon=dencon, overden=overden)
		
		self.centers_indices = np.asarray(self.centers_grid >= self.vote_threshold).nonzero()
		self.C_weights = self.centers_grid[self.centers_indices]
		if self.printout:
			precut = len(self.centers_grid[self.centers_grid!=0])
			postcut = len(self.C_weights)
			print('Number of found centers before vote cut:', precut)
			print('Number of found centers after vote cut:', postcut)
		delattr(self, 'centers_grid')

		# calculates center coords to be exactly at the center of their respective bins
		centers_bin_coords = np.array([(self.density_grid_edges[i][:-1] + self.density_grid_edges[i][1:]) / 2 
			for i in range(len(self.density_grid_edges))])
		delattr(self, 'density_grid_edges')
		C_xyzs = np.array([centers_bin_coords[i][self.centers_indices[i]] for i in range(len(self.centers_indices))])
		self.C_ra, self.C_dec, self.C_redshift, _ = cartesian2sky(*C_xyzs, self.LUT_redshifts, 
																self.G_ra.min(), self.G_ra.max())
		
		# outputs centers catalog in skycoords+weights to out folder
		savename = self.savename + f'found_centers_r_{self.kernel_radius}_cut_{self.vote_threshold}.fits'
		save_data_weighted(savename, self.C_ra, self.C_dec, self.C_redshift, self.C_weights)



	"""
	DEV IN PROGRESS
	"""
	def refine_centers(self, refinement_args):
		"""
		Refines the locations of found centers.
		"""

		kernel_type = refinement_args[0]
		kernel_args = [float(refinement_args[1])]
		# grid_spacing = refinement_args[2]

		# makes the kernel for scanning over the density grid
		# kernel_grid = self._kernel()
		kernel = Kernel(kernel_type, self.kernel_radius, self.grid_spacing,
			self.printout, self.show_kernel, *kernel_args)

		# for refinement process, the following allows kernel jittering
		# box around each center is 10% larger than the kernel outer radius
		self.centerbox_r_upper_bound_idx_units = int(np.round(kernel.\
			kernel_r_idx_units_upper_bound*1.1))


		# print(self.centers_indices)

		for idxs in zip(*self.centers_indices):

			# note down original indexes because will have to do a coord
			# translation in the end after the convolution

			# print(idxs)

			# need the box of galaxies around this center
			# these calculations deal with edge cases
			(xlo,xhi),(ylo,yhi),(zlo,zhi) = tuple(
			(max(0, idxs[i]-self.centerbox_r_upper_bound_idx_units), 
			min(idxs[i]+self.centerbox_r_upper_bound_idx_units, self.density_grid.shape[i]))
			for i in range(len(idxs)))

			# print(centerbox_bounds_idx_units)
			# print(xlo,xhi,ylo,yhi,zlo,zhi)

			densitybox = self.density_grid[xlo:xhi,ylo:yhi,zlo:zhi]
			# print(densitybox.shape)

			centerbox = np.round(fftconvolve(densitybox, kernel.get_grid(), mode='same'))

			# print(centerbox.shape)


			# break








