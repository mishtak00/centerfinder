import os
from argparse import ArgumentParser
from centerfinder import CenterFinder



def main():
	parser = ArgumentParser(description="~~~~~~~~~~~~~~~~~ ( * ) Center Finder ( * ) ~~~~~~~~~~~~~~~~~")
	parser.add_argument('file', metavar='GALAXY_FILE', type=str, 
		help='Name of fits file with the galaxy data.')
	parser.add_argument('-r', '--kernel_radius', type=float, default=108., 
		help='If this argument is present, the kernel radius will be set to the value entered as argument.')
	parser.add_argument('-t', '--vote_threshold', type=float, default=0.,
		help='If this argument is present, centers with number of votes smaller than given argument will be discarded from .fits output.')
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

	cf = CenterFinder(args.file, args.kernel_radius, args.vote_threshold, params_file=args.params_file, save=args.save)
	cf.find_centers()


if __name__ == '__main__':
	main()