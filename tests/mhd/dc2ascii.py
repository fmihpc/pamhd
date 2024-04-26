#! /usr/bin/env python3
'''
Converts output of MHD test program of PAMHD to ASCII format.

Copyright 2015, 2016 Ilja Honkonen
Copyright 2024 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Author(s): Ilja Honkonen
'''

from importlib.util import module_from_spec, spec_from_file_location
from os.path import dirname, isfile, join
from sys import modules
spec = spec_from_file_location('common', join(dirname(__file__), 'common.py'))
common = module_from_spec(spec)
modules['common'] = common
spec.loader.exec_module(common)

try:
	import numpy
except Exception as e:
	print("Couldn't import numpy, try pip install --user numpy: ", e)
	exit(1)


def convert(inname, verbose):
	if not isfile(inname):
		raise Exception('Given path (' + inname + ') is not a file.')

	infile = open(inname, 'rb')

	if verbose:
		print('Converting file', inname)

	metadata = common.get_metadata(infile)
	metadata['cells'] = common.get_cells(infile, metadata['total_cells'])
	sim_data_ = common.get_cell_data(infile, metadata, range(len(metadata['cells'])))
	sim_data = dict()
	for i in range(len(metadata['cells'])):
		sim_data[metadata['cells'][i]] = sim_data_[i]

	if metadata['file_version'] != 3:
		exit('Unsupported file version: ' + str(metadata['file_version']))
	if verbose:
		print('Simulation step:', metadata['sim_step'])
		print('Simulation time:', metadata['sim_time'])

	if metadata['geometry_id'] != 1:
		raise Exception('Unsupported geometry')

	if verbose:
		print('Number of simulation cells:', metadata['total_cells'])

	outname = inname.replace('.dc', '.txt')
	outfile = open(outname, 'w')

	if verbose:
		print('Number of variables:', len(metadata['var_data_start']), '\nVariables:')
		for name in metadata['var_data_start']:
			print(name, 'start in file:', metadata['var_data_start'])

	outfile.write(
		  "# MHD data created by PAMHD\n"
		+ "# Adiabatic index: " + str(metadata['adiabatic_index'])
		+ "\n# Vacuum permeability: " + str(metadata['vacuum_permeability'])
		+ "\n# Simulation time: " + str(metadata['sim_time'])
		+ "\n# Simulation step: " + str(metadata['sim_step']) + "\n"
		+ "# Format of each line (lines are ordered randomly):\n"
		+ "# x, y and z coordinates of cell's center (m)\n")
	if 'mhd' in metadata['var_data_start']: outfile.write(
		"# number density (#/m^3)\n"
		+ "# x, y and z components of velocity (m/s)\n"
		+ "# thermal pressure (Pa)\n"
		+ "# x, y and z components of volume magnetic field (T)\n")
	if 'primary' in metadata['var_data_start']: outfile.write(
		"# primary face info (-x,+x,-y...,+z)\n"
		+ "# primary edge info (x,-y,-z; x,-y,+z; ...; z,+x,+y)\n")
	if 'bgB' in metadata['var_data_start']: outfile.write("# background magnetic field on cell faces (-x,+x,-y,...,+z) with x,y,z components on each face (T)\n")
	if 'divfaceB' in metadata['var_data_start']: outfile.write("# divergence of face magnetic field (T/m)\n")
	if 'edgeE' in metadata['var_data_start']: outfile.write("# edge-directed electric field on cell edges (V/m, x,-y,-z; x,-y,+z; ...; z,+x,+y)\n")
	if 'faceB' in metadata['var_data_start']: outfile.write("# normal component of magnetic field on faces (T, -x,+x,-y,...,+z)\n")
	if 'mpi' in metadata['var_data_start']: outfile.write("# MPI rank of cell owner\n")
	if 'mhd info' in metadata['var_data_start']: outfile.write("# MHD solver information\n")
	if 'ref lvls' in metadata['var_data_start']: outfile.write("# Minimum and maximum target refinement level\n")

	for i in range(len(metadata['cells'])):
		cell_id = metadata['cells'][i]
		cell_center = common.get_cell_center(metadata, cell_id)

		outfile.write(
			str(cell_center[0]) + ' '
			+ str(cell_center[1]) + ' '
			+ str(cell_center[2]) + ' ')

		if 'mhd' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['mhd'] + i*8*8, 0)
			mas, mom, nrj, mag = numpy.fromfile(
				infile,
				dtype = 'double, 3double, double, 3double',
				count = 1
			)[0]

			kin_nrj = (mom[0]**2 + mom[1]**2 + mom[2]**2) / 2 / mas
			mag_nrj = (mag[0]**2 + mag[1]**2 + mag[2]**2) / 2 / metadata['vacuum_permeability']
			pressure = (nrj - kin_nrj - mag_nrj) * (metadata['adiabatic_index'] - 1)
			outfile.write(
				str(mas / metadata['proton_mass']) + ' '
				+ str(mom[0]/mas) + ' '
				+ str(mom[1]/mas) + ' '
				+ str(mom[2]/mas) + ' '
				+ str(pressure) + ' '
				+ str(mag[0]) + ' '
				+ str(mag[1]) + ' '
				+ str(mag[2]) + ' ')

		if 'primary' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['primary'] + i*(6+12), 0)
			p = numpy.fromfile(infile, dtype = '18u1', count = 1)[0]
			for i in range(len(p)):
				outfile.write(str(p[i]) + ' ')

		if 'bgB' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['bgB'] + i*3*6*8, 0)
			B = numpy.fromfile(infile, dtype = '18double', count = 1)[0]
			for i in range(len(B)):
				outfile.write(str(B[i]) + ' ')

		if 'divfaceB' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['divfaceB'] + i*8, 0)
			B = numpy.fromfile(infile, dtype = 'double', count = 1)[0]
			outfile.write(str(B) + ' ')

		if 'edgeE' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['edgeE'] + i*12*8, 0)
			E = numpy.fromfile(infile, dtype = '12double', count = 1)[0]
			for i in range(len(E)):
				outfile.write(str(E[i]) + ' ')

		if 'faceB' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['faceB'] + i*6*8, 0)
			B = numpy.fromfile(infile, dtype = '6double', count = 1)[0]
			for i in range(len(B)):
				outfile.write(str(B[i]) + ' ')

		if 'rank' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['rank'] + i*4, 0)
			rank = numpy.fromfile(infile, dtype = 'intc', count = 1)[0]
			outfile.write(str(rank) + ' ')

		if 'mhd info' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['mhd info'] + i*4, 0)
			s = numpy.fromfile(infile, dtype = 'uintc', count = 1)[0]
			outfile.write(str(s) + ' ')

		if 'ref lvls' in metadata['var_data_start']:
			infile.seek(metadata['var_data_start']['ref lvls'] + i*8, 0)
			r = numpy.fromfile(infile, dtype = '2intc', count = 1)[0]
			outfile.write(str(r[0]) + ' ' + str(r[1]) + ' ')

		outfile.write('\n')


if __name__ == '__main__':
	import sys
	verbose = False
	for arg in sys.argv[1:]:
		if arg == '--verbose':
			verbose = True
		else:
			convert(arg, verbose)
