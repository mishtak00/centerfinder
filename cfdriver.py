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
from centerfinder import CenterFinder



def main():
	parser = ArgumentParser(description="~~~~~~~~~~~~~~~~~ ( * ) Center Finder ( * ) ~~~~~~~~~~~~~~~~~")
	parser.add_argument('file', metavar='GALAXY_FILE', type=str, 
		help='Name of fits file with the galaxy data.')
	parser.add_argument('-r', '--kernel_radius', type=float, 
		help='If this argument is present, the kernel radius will be set to the value entered as argument.')
	parser.add_argument('-t', '--vote_threshold', type=float,
		help='If this argument is present, centers with number of votes smaller than given argument will be discarded from .fits output.')
	parser.add_argument('-w', '--weighted_input', action='store_true',
		help='If this argument is present, CenterFinder will try to read a fourth column from input data and interpret said values as weights.')
	back_or_over = parser.add_mutually_exclusive_group()
	back_or_over.add_argument('-c', '--density_contrast', action='store_true',
		help='If this argument is present, the CenterFinder will subtract the background from the galaxy density grid before voting.')
	back_or_over.add_argument('-o', '--overdensity', action='store_true',
		help='If this argument is present, the CenterFinder will subtract average density from the galaxy density grid before voting.')
	parser.add_argument('-p', '--params_file', type=str, default='params.json', 
		help='If this argument is present, the cosmological parameters will be loaded from file given as argument.')
	parser.add_argument('-s', '--save', action='store_true', 
		help='If this argument is present, grids and .fits output will be automatically saved to an \'out\' folder.')
	parser.add_argument('-v', '--verbose', action='store_true', 
		help='If this argument is present, the progress of center-finder will be printed out to standard output.')
	args = parser.parse_args()

	# deletes the .fits extension
	# allows for other '.'s in the args.file string
	filename = '.'.join(args.file.split('.')[:-1])

	if args.save:
		try:
			os.mkdir('out_{}'.format(filename))
		except FileExistsError:
			pass

	cf = CenterFinder(args.file, args.weighted_input, 
		args.params_file, args.save, args.verbose)
	if args.kernel_radius is not None:
		cf.set_kernel_radius(args.kernel_radius)
	if args.vote_threshold is not None:
		cf.set_vote_threshold(args.vote_threshold)

	cf.find_centers(dencon=args.density_contrast,
					overden=args.overdensity)



if __name__ == '__main__':
	main()
