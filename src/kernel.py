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

		if self.printout:
			print('Constructing kernel...')

		if self.type=='step':
			# interprets args[0] as step function thickness in h^-1Mpc
			self.thickness = self.args[0] if len(self.args)>0 else 1
			self.grid = self._step(self.thickness)
		# elif self.type=='gaussian':


		if self.printout:
			print('Kernel constructed successfully...')
			print('Number of kernel bins that contain surface:', len(self.grid[self.grid!=0]))
			print('Number of empty kernel bins:', len(self.grid[self.grid==0]))

		# this is here for sanity checks
		# shows the kernel in 3D with blue disks in nonzero kernel bins
		if self.plot:
			import matplotlib.pyplot as plt
			from mpl_toolkits.mplot3d import Axes3D
			color = 'cornflowerblue'
			fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
			ax.scatter(*np.where(self.grid!=0), c=color)
			plt.show()



	def get_grid(self):
		return self.grid


	def _step(self, thickness: float) -> np.ndarray:
		# calculates the kernel inscribed radius in index units.
		# for a kernel with radius 100 and grid_spacing 10, the radius
		# is 10 in idx unitx. this puts the value of the function at
		# r = 100 in the 11th bin of the 1D function.
		kernel_r_idx_units = int(np.ceil(self.radius / self.grid_spacing))
		kernel_r_idx_units_upper_bound = kernel_r_idx_units + 0.5 * thickness
		kernel_r_idx_units_lower_bound = kernel_r_idx_units - 0.5 * thickness

		# calculates the number of bins in each dimensional axis
		# this calculation ensures an odd numbered gridding so that
		# the kernel construction has a distinct central bin on any given run
		kernel_bin_count = 2 * kernel_r_idx_units + 1

		# central bin index, since the kernel is a cube this can just be one int
		kernel_center_idx = int(kernel_bin_count / 2)
		kernel_center = np.array([kernel_center_idx, ] * 3)

		# evaluates 1D function that's gonna be rotated below
		# defined over the entire radius of the circumscribed sphere of the kernel cube
		circumscribed_r_idx_units = int(np.ceil(3**.5 * kernel_r_idx_units))
		

		if self.printout:
			print('Kernel radius in index units:', kernel_r_idx_units)
			print('Kernel radius upper bound:', kernel_r_idx_units_upper_bound)
			print('Kernel radius lower bound:', kernel_r_idx_units_lower_bound)
			print('Kernel bin count:', kernel_bin_count)

		step_func = np.array([
			1 if i>=kernel_r_idx_units_lower_bound 
			and i<kernel_r_idx_units_upper_bound 
			else 0
			for i in range(circumscribed_r_idx_units)])

		# normalization = 1 / integral
		integral = sum(step_func)
		step_func = step_func / integral

		# print(step_func)

		return np.array([[[
			step_func[int(np.linalg.norm(np.array([i,j,k])-kernel_center))]
			for k in range(kernel_bin_count)]
			for j in range(kernel_bin_count)]
			for i in range(kernel_bin_count)])





# if __name__ == '__main__':
# 	Kernel('step', 110, 5, True, True, 2)


