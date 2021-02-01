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
from argparse import ArgumentParser
from src.centerfinder import CenterFinder



def main():

	print('\n\n\nPlease reference original publication arXiv:2008.12793 '\
		'when using this software for publishing/redistribution.\n\n\n')
	
	parser = ArgumentParser(description=
		'~~~~~~~~~~~~~~~~~ ( * ) Center Finder ( * ) ~~~~~~~~~~~~~~~~~')

	# these read input and parameter files
	parser.add_argument('file', metavar='GALAXY_FILE', type=str, 
		help='Name of .fits file with the input data.')
	parser.add_argument('-p', '--params_file', type=str, default='params.json', 
		help='Sets custom hyperparameters file.')
	
	# these define kernel behavior
	parser.add_argument('-r', '--kernel_radius', type=float, help='Sets kernel radius.')
	parser.add_argument('--show_kernel', action='store_true', help='Shows 1D kernel plot.')
	kernel_types = parser.add_mutually_exclusive_group()
	kernel_types.add_argument('-e', '--step_kernel', nargs='*', type=float,
		help='Fits a step function to the kernel at kernel radius.')
	kernel_types.add_argument('-g', '--gaussian_kernel', nargs=1, type=float,
		help='Fits a gaussian function to the kernel at kernel radius.')
	kernel_types.add_argument('-a', '--wavelet_kernel', nargs=1, type=float, 
		help='Fits a wavelet function to the kernel at kernel radius.')
	kernel_types.add_argument('-u', '--custom_kernel', nargs=1, type=str, 
		help='Fits given custom array to kernel radially.')

	# these define behavior of density grid
	parser.add_argument('-t', '--vote_threshold', type=float, 
		help='Centers with number of votes smaller than given argument '\
		'will be discarded from .fits output.')
	parser.add_argument('-w', '--weighted_input', action='store_true',
		help='CenterFinder will try to read a fourth column from input data '\
		'and interpret said values as weights.')
	con_or_over = parser.add_mutually_exclusive_group()
	con_or_over.add_argument('-c', '--density_contrast', nargs='*',
		help='CenterFinder will subtract the background from the galaxy '\
		'density grid before voting. It will set negative weights to 0 if '\
		'anything is entered after -c.')
	con_or_over.add_argument('-o', '--overdensity', action='store_true',
		help='CenterFinder will subtract average density from the galaxy '\
		'density grid before voting.')
	
	# ancillary behaviors
	parser.add_argument('-s', '--save', action='store_true', 
		help='Grids and .fits output will be automatically saved to an \'out\' folder.')
	parser.add_argument('-v', '--verbose', action='store_true', 
		help='The progress of CenterFinder will be printed out to standard output.')
	parser.add_argument('-l', '--plot_slice', nargs='*', type=float,
		help='CenterFinder will plot a slice of the centers found with '\
		'galaxies that voted for them superimposed.')
	parser.add_argument('--plot_coord_hist', type=str,
		help='CenterFinder will plot a histogram of the selected coordinate for '\
		'both galaxies and centers. Choices: RA, DEC, Z, R')

	# refinement process args
	parser.add_argument('--refine', nargs='*')

	args = parser.parse_args()

	# deletes the .fits extension and
	# allows for other '.'s in the args.file string
	filename = '.'.join(args.file.split('.')[:-1])
	try:
		os.mkdir('out_{}'.format(filename))
	except FileExistsError:
		pass

	# creates and customizes instance of CenterFinder object
	cf = CenterFinder(args.file, args.weighted_input, 
		args.params_file, args.save, args.verbose)
	if args.kernel_radius is not None:
		cf.set_kernel_radius(args.kernel_radius)
	if args.show_kernel:
		cf.set_show_kernel(args.show_kernel)
	if args.step_kernel is not None:
		cf.set_kernel_type('step', args.step_kernel)
	elif args.gaussian_kernel is not None:
		cf.set_kernel_type('gaussian', args.gaussian_kernel)
	elif args.wavelet_kernel is not None:
		cf.set_kernel_type('wavelet', args.wavelet_kernel)
	elif args.custom_kernel is not None:
		cf.set_kernel_type('custom', args.custom_kernel)
	if args.vote_threshold is not None:
		cf.set_vote_threshold(args.vote_threshold)

	# runs the centerfinding algorithm
	if args.density_contrast is not None:
		do_dencon = True
		if len(args.density_contrast)==0:
			keep_neg_wts = True
		else:
			keep_neg_wts = False
	else:
		do_dencon = False
		keep_neg_wts = False
	dencon_args = (do_dencon, keep_neg_wts)
	cf.find_centers(dencon=dencon_args,
					overden=args.overdensity)

	# plotting options
	if args.plot_slice is not None:
		from src.plotter import _plot_slice
		_plot_slice(cf, args.plot_slice)
	if args.plot_coord_hist is not None:
		from src.plotter import _plot_coord_hist
		_plot_coord_hist(cf, args.plot_coord_hist)


	# jittering DEV
	# if args.refine is not None:
	# 	cf.refine_centers(args.refine)


if __name__ == '__main__':
	main()




