from math import pi, exp
import numpy as np


class Kernel:


	def __init__(self, type: str, radius: float, grid_spacing: float, 
		printout: bool, plot: bool, *args):

		self.type = type
		self.radius = radius
		self.args = args
		self.grid_spacing = grid_spacing
		self.plot = plot
		self.printout = printout

		if self.type=='step':
			# interprets args[0] as step function thickness in h^-1Mpc
			self.thickness = self.args[0] if len(self.args)>0 else 1
		elif self.type=='gaussian':
			# interprets args[0] as stdev of gaussian with mean at kernel rad
			self.stdev = self.args[0]
		elif self.type=='wavelet':
			# interprets args[0] as width of wavelet with scale at kernel rad
			self.width = self.args[0]
		elif self.type=='custom':
			# interprets args[0] as the file containing the custom array
			self.source = self.args[0]

		if self.printout:
			print(f'Constructing {self.type} kernel...')

		self.grid, self.kernel_r_idx_units_upper_bound = self._make_grid()

		if self.printout:
			print('Kernel constructed successfully...')
			print('Number of nonzero kernel bins:', len(self.grid[self.grid!=0]))
			print('Number of empty kernel bins:', len(self.grid[self.grid==0]))

		# # this is here for sanity checks
		# # shows the kernel in 3D with blue disks in nonzero kernel bins
		# if self.plot:
		# 	import matplotlib.pyplot as plt
		# 	from mpl_toolkits.mplot3d import Axes3D
		# 	fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
		# 	ax.scatter(*np.where(self.grid!=0), c='cornflowerblue')
		# 	plt.show()



	def get_grid(self):
		return self.grid



	def _calculate_settings(self):
		# calculates the kernel inscribed radius in index units.
		# for a kernel with radius 100 and grid_spacing 10, the radius
		# is 10 in idx unitx. this puts the value of the function at
		# r = 100 in the 11th bin of the 1D function.
		kernel_r_idx_units = int(np.ceil(self.radius / self.grid_spacing))
		# calculates circumscribed sphere radius for easy eval of 1D func
		circumscribed_r_idx_units = int(np.ceil(3**.5 * kernel_r_idx_units))
		# calculates the number of bins in each dimensional axis
		# this calculation ensures an odd numbered gridding so that
		# the kernel construction has a distinct central bin on any given run
		kernel_bin_count = 2 * kernel_r_idx_units + 1
		# central bin index, since the kernel is a cube this can just be one int
		kernel_center_idx = int(kernel_bin_count / 2)
		kernel_center = np.array([kernel_center_idx, ] * 3)

		return (
			kernel_r_idx_units, circumscribed_r_idx_units,
			kernel_bin_count, kernel_center_idx, kernel_center)


	def _calculate_settings_custom(self, func_arr: np.array):
		# look at _calculate_settings() for descriptions of below calculations
		kernel_r_idx_units = len(func_arr)
		circumscribed_r_idx_units = int(np.ceil(3**.5 * kernel_r_idx_units))
		kernel_bin_count = 2 * kernel_r_idx_units + 1
		kernel_center_idx = int(kernel_bin_count / 2)
		kernel_center = np.array([kernel_center_idx, ] * 3)
		return (
			kernel_r_idx_units, circumscribed_r_idx_units,
			kernel_bin_count, kernel_center_idx, kernel_center)



	def _make_step(self, kernel_r_idx_units: int, circumscribed_r_idx_units: int):
		# transforms kernel thickness to index units
		thickness_idx_units = int(np.ceil(self.thickness / self.grid_spacing))
		kernel_r_idx_units_upper_bound = kernel_r_idx_units + 0.5 * thickness_idx_units
		kernel_r_idx_units_lower_bound = kernel_r_idx_units - 0.5 * thickness_idx_units

		step_func = np.array([
			1 if i>=kernel_r_idx_units_lower_bound 
			and i<kernel_r_idx_units_upper_bound 
			else 0
			for i in range(circumscribed_r_idx_units)])

		# normalization = 1 / integral
		step_func = step_func / np.sum(step_func) # normalizes

		return step_func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound



	def _make_gaussian(self, kernel_r_idx_units: int, circumscribed_r_idx_units: int):
		# transforms stdev to index units
		stdev_idx_units = int(np.ceil(self.stdev / self.grid_spacing))
		# 99.7% of all data in a normal dist is within 3 std devs of mean
		kernel_r_idx_units_upper_bound = kernel_r_idx_units + 3*stdev_idx_units
		kernel_r_idx_units_lower_bound = kernel_r_idx_units - 3*stdev_idx_units

		# calculates 1/(stdev root(2pi)) e^(- (x-mean)^2 / (2 stdev^2))
		over_sroot2pi = 1. / (stdev_idx_units * (2.*pi)**.5)
		minus_over_2ssquared = -1. / (2. * stdev_idx_units**2)
		gauss_func = np.array([
			(lambda mean,stdev,x: 
				over_sroot2pi * exp( minus_over_2ssquared * (x-mean)**2 ))
			(kernel_r_idx_units, stdev_idx_units,i)
			if i>=kernel_r_idx_units_lower_bound 
			and i<kernel_r_idx_units_upper_bound 
			else 0
			for i in range(circumscribed_r_idx_units)])
		gauss_func = gauss_func / np.sum(gauss_func) # normalizes

		return gauss_func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound



	def _make_wavelet(self, kernel_r_idx_units: int, circumscribed_r_idx_units: int):
		# transforms wavelet width parameter to index units
		width_idx_units = int(np.ceil(self.width / self.grid_spacing))
		kernel_r_idx_units_upper_bound = kernel_r_idx_units + 2*width_idx_units
		kernel_r_idx_units_lower_bound = kernel_r_idx_units - 2*width_idx_units

		# calculates psi(x) = 1/(4 pi x^2) ( 2 B_3(2(x-R)/s) - B_3((x-R)/s) )
		# where B_3(x) = 1/12 (|x-2|^3 - 4|x-1|^3 + 6|x|^3 - 4|x+1|^3 + |x+2|^3)
		# R is scale = kernel radius and s is width, both are in index units
		@np.vectorize
		def _B_3(x):
			return (abs(x-2)**3 - 4*abs(x-1)**3 + 6*abs(x)**3 - 4*abs(x+1)**3 + abs(x+2)**3) / 12
		@np.vectorize
		def _psi(x):
			y = (x - kernel_r_idx_units) / width_idx_units
			return (2*_B_3(2*y) - _B_3(y)) / (4*pi*x**2)
		# wave_func = np.array([psi(i) for i in range(circumscribed_r_idx_units)])
		# wave_func = _psi(np.arange(circumscribed_r_idx_units))
		wave_func = np.array([
			_psi(i)
			if i>=kernel_r_idx_units_lower_bound 
			and i<kernel_r_idx_units_upper_bound
			else 0
			for i in range(circumscribed_r_idx_units)])
		wave_func = wave_func / np.sum(abs(wave_func)) # normalizes

		return wave_func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound



	def _make_custom(self):
		# reads custom array from file passed as argument to Kernel
		custom_func = np.load(self.source)
		# kernel radius here is the whole domain of the user defined function
		# multiplied by the grid_spacing, so a function defined over 22 indices with
		# a grid_spacing 5Mpc/h infers that the kernel_radius is 110Mpc/h
		kernel_r_idx_units_lower_bound = 0
		kernel_r_idx_units_upper_bound = len(custom_func)

		return custom_func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound



	def _make_grid(self) -> np.ndarray:

		if self.type=='step' or self.type=='gaussian' or self.type=='wavelet':
			kernel_r_idx_units, circumscribed_r_idx_units,\
			kernel_bin_count, kernel_center_idx, kernel_center = self._calculate_settings()
			

		# evaluates 1D function that's gonna be rotated below
		# defined over the entire radius of the circumscribed sphere of the kernel cube
		if self.type=='step':
			func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound = \
				self._make_step(kernel_r_idx_units, circumscribed_r_idx_units)
		elif self.type=='gaussian':
			func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound = \
				self._make_gaussian(kernel_r_idx_units, circumscribed_r_idx_units)
		elif self.type=='wavelet':
			func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound = \
				self._make_wavelet(kernel_r_idx_units, circumscribed_r_idx_units)
		elif self.type=='custom':
			# the upper bound radius is just the radius here
			func, kernel_r_idx_units_upper_bound, kernel_r_idx_units_lower_bound = self._make_custom()
			kernel_r_idx_units, circumscribed_r_idx_units,\
			kernel_bin_count, kernel_center_idx, kernel_center = self._calculate_settings_custom(func)
			# pads function with 0s from the user-fed radius to the circurmscribed radius
			func = np.array([func[i] if i<kernel_r_idx_units else 0 for i in range(circumscribed_r_idx_units)])

		# sets

		if self.printout:
			print('Kernel radius in index units:', kernel_r_idx_units)
			print('Kernel radius upper bound:', kernel_r_idx_units_upper_bound)
			print('Kernel radius lower bound:', kernel_r_idx_units_lower_bound)
			print('Kernel bin count (side length of grid in idx units):', kernel_bin_count)
			print('Integral of 1D kernel function:', np.sum(abs(func)))

		if self.plot:
			self._plot_1d_func(func, circumscribed_r_idx_units)

		return Kernel._sphericize(func, kernel_center, kernel_bin_count), kernel_r_idx_units_upper_bound



	@staticmethod
	def _sphericize(func: np.array, center: np.array, bin_count: int):
		return np.array([[[
			func[int(np.linalg.norm(np.array([i,j,k])-center))]
			for k in range(bin_count)]
			for j in range(bin_count)]
			for i in range(bin_count)])



	def _plot_1d_func(self, func: np.array, outer_bound: float):
		import matplotlib.pyplot as plt
		bins = self.grid_spacing * np.arange(0.,outer_bound)
		plt.bar(bins, func, width=self.grid_spacing, edgecolor='black')
		plt.title(f'Bar plot of discrete kernel function in 1D: {self.type.upper()}')
		plt.xlabel('Radial distance from kernel center $[h^{-1}Mpc]$')
		plt.show()







