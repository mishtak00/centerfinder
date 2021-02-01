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

import json
import numpy as np
from astropy.io import fits
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline



def load_hyperparameters(params_file: str):
	"""Returns hyperparameters loaded from json file. """
	with open(params_file, 'r') as params:
		hp = json.load(params)
		c = hp['c']  # km/s
		H0 = hp['H0']
		c_over_H0 = c / H0
		h_ = H0 / 100.
		Omega_M = hp['Omega_M']
		Omega_K = hp['Omega_K']
		if 'Omega_L' not in hp:
			Omega_L = 1 - Omega_M - Omega_K
		else:
			Omega_L = hp['Omega_L']
		grid_spacing = hp['grid_spacing']  # h-1Mpc
	cosmology = h_, c_over_H0, Omega_M, Omega_K, Omega_L
	return cosmology, grid_spacing


def load_data(filename: str) -> (np.recarray,)*3:
	with fits.open(filename) as catalog_data:
		catalog_data = catalog_data[1].data
		ra = catalog_data['ra']
		dec = catalog_data['dec']
		redshift = catalog_data['z']
	return ra, dec, redshift


def load_data_weighted(filename: str) -> (np.recarray,)*4:
	with fits.open(filename) as catalog_data:
		catalog_data = catalog_data[1].data
		ra = catalog_data['ra']
		dec = catalog_data['dec']
		redshift = catalog_data['z']
		weight = catalog_data['wts']
	return ra, dec, redshift, weight


def save_data(filename: str, ra: np.array, dec: np.array, z: np.array):
	racol = fits.Column(name='ra', array=ra, format='E')
	deccol = fits.Column(name='dec', array=dec, format='E')
	zcol = fits.Column(name='z', array=z, format='E')
	t = fits.BinTableHDU.from_columns([racol, deccol, zcol])
	t.writeto(filename, overwrite=True)


def save_data_weighted(filename: str, ra: np.array, dec: np.array, z: np.array, w: np.array):
	# 'E' is the fits format for single-precision floating point
	racol = fits.Column(name='ra', array=ra, format='E')
	deccol = fits.Column(name='dec', array=dec, format='E')
	zcol = fits.Column(name='z', array=z, format='E')
	wcol = fits.Column(name='wts', array=w, format='E')
	t = fits.BinTableHDU.from_columns([racol, deccol, zcol, wcol])
	t.writeto(filename, overwrite=True)


def z2r(z: np.array, cosmology: tuple) -> float:
	"""Transforms observed redshift to radial distance. """
	h_, c_over_H0, Omega_M, Omega_K, Omega_L = cosmology
	# multiply by lil h to transform from Mpc to h-1Mpc
	const = h_ * c_over_H0
	return const * integrate.quad(lambda z_: 
		(Omega_M*(1+z_)**3 + Omega_K*(1+z_)**2 + Omega_L)**-0.5, 0, z)[0]


def interpolate_r_z(redshift_min: float, redshift_max: float, cosmology: tuple):
	# TODO: automate getting the number of spacings instead of setting it to 100
	# Creates tick values for the given range of redshifts
	redshift_ticks = np.linspace(redshift_min, redshift_max, 100)
	radii_ticks = np.array([z2r(z, cosmology) for z in redshift_ticks])
	# creates lookup tables with interpolated radii values
	LUT_radii = InterpolatedUnivariateSpline(redshift_ticks, radii_ticks)
	LUT_redshifts = InterpolatedUnivariateSpline(radii_ticks, redshift_ticks)
	return LUT_radii, LUT_redshifts


def sky2cartesian(ra: np.array, dec: np.array, redshift: np.array, LUT_radii) -> (np.array, np.array, np.array):
	# note that polar angle = pi/2 - declination
	# sin(polar) = cos(90-polar) = cos(dec)
	# cos(polar) = sin(90-polar) = sin(dec)
	radii = np.array([LUT_radii(z) for z in redshift])
	xs = radii * np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(ra))
	ys = radii * np.cos(np.deg2rad(dec)) * np.sin(np.deg2rad(ra))
	zs = radii * np.sin(np.deg2rad(dec))
	return xs, ys, zs


def cartesian2sky(xs: np.array, ys: np.array, zs: np.array, LUT_redshifts, 
	ramin: float, ramax: float) -> (np.array,)*4:
	radii = (xs ** 2 + ys ** 2 + zs ** 2) ** 0.5
	redshift = np.array(LUT_redshifts(radii))
	dec = np.rad2deg(np.arcsin(zs / radii))
	# TODO: fix this hack
	rafix = np.round((ramin+ramax)/2 //90) *90
	ra = np.rad2deg(np.arctan(ys / xs)) + rafix
	return ra, dec, redshift, radii


