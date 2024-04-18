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

import os

try:
	import numpy
except Exception as e:
	print("Couldn't import numpy, try pip install --user numpy: ", e)
	exit(1)


def convert(inname, verbose):
	if not os.path.isfile(inname):
		raise Exception('Given path (' + inname + ') is not a file.')

	infile = open(inname, 'rb')

	if verbose:
		print('Converting file', inname)

	# read simulation header, given by source/mhd/save.hpp
	file_version = numpy.fromfile(infile, dtype = 'uint64', count = 1)
	if file_version != 3:
		exit('Unsupported file version: ' + str(file_version))

	sim_step = numpy.fromfile(infile, dtype = 'uint64', count = 1)[0]
	if verbose:
		print('Simulation step:', sim_step)

	sim_time, adiabatic_index, proton_mass, vacuum_permeability \
		= numpy.fromfile(infile, dtype = '4double', count = 1)[0]
	if verbose:
		print('Simulation time:', sim_time)

	nr_offsets = numpy.fromfile(infile, dtype = 'u1', count = 1)[0]
	var_offsets = numpy.fromfile(infile, dtype = 'uint64', count = nr_offsets)

	# check file endiannes
	endianness = numpy.fromfile(infile, dtype = 'uint64', count = 1)[0]
	if endianness != numpy.uint64(0x1234567890abcdef):
		raise Exception(
			'Wrong endiannes in given file, expected ' + str(hex(0x1234567890abcdef)) \
			+ ' but got ' + str(hex(endianness))
		)

	# number of refinement level 0 cells in each dimension
	ref_lvl_0_cells = numpy.fromfile(infile, dtype = '3uint64', count = 1)[0]
	#print(ref_lvl_0_cells)

	# maximum refinement level of grid cells
	max_ref_lvl = numpy.fromfile(infile, dtype = 'intc', count = 1)[0]
	if max_ref_lvl > numpy.uint32(0):
		raise Exception('Refinement level > 0 not supported')

	# length of every cells' neighborhood in cells of identical size
	neighborhood_length = numpy.fromfile(infile, dtype = 'uintc', count = 1)[0]
	#print(neighborhood_length)

	# whether grid is periodic each dimension (0 == no, 1 == yes)
	periodicity = numpy.fromfile(infile, dtype = '3uint8', count = 1)[0]
	#print(periodicity)

	geometry_id = numpy.fromfile(infile, dtype = 'intc', count = 1)[0]
	if geometry_id != numpy.int32(1):
		raise Exception('Unsupported geometry')

	# starting coordinate of grid
	grid_start = numpy.fromfile(infile, dtype = '3double', count = 1)[0]
	#print(grid_start)

	# length of cells of refinement level 0
	lvl_0_cell_length = numpy.fromfile(infile, dtype = '3double', count = 1)[0]
	#print(lvl_0_cell_length)

	# total number of cells in grid
	total_cells = numpy.fromfile(infile, dtype = 'uint64', count = 1)[0]
	if verbose:
		print('Number of simulation cells:', total_cells)

	# id of each cell
	cell_ids = numpy.fromfile(infile, dtype = 'uint64', count = total_cells)

	outname = inname.replace('.dc', '.txt')
	outfile = open(outname, 'w')
	data = dict()
	name_len = 8
	for offset in var_offsets:
		infile.seek(offset, 0)
		name = str(numpy.fromfile(infile, dtype = 'byte', count = name_len), 'ascii').strip()
		data[name] = int(offset + name_len)

	if verbose:
		print('Number of variables:', len(data.keys()), '\nVariables:')
		for name in data:
			print(name, 'start in file:', data[name])

	outfile.write(
		  "# MHD data created by PAMHD\n"
		+ "# All vectors in Cartesian GSE coordinate system\n"
		+ "# Format of each line (lines are ordered randomly):\n"
		+ "# x, y and z coordinates of cell's center (m)\n")
	if 'mhd' in data: outfile.write(
		"# number density (#/m^3)\n"
		+ "# x, y and z components of velocity (m/s)\n"
		+ "# thermal pressure (Pa)\n"
		+ "# x, y and z components of volume magnetic field (T)\n")
	if 'primary' in data: outfile.write(
		"# primary face info (-x,+x,-y...,+z)\n"
		+ "# primary edge info (x,-y,-z; x,-y,+z; ...; z,+x,+y)\n")
	if 'bgB' in data: outfile.write("# background magnetic field on cell faces (-x,+x,-y,...,+z) with x,y,z components on each face (T)\n")
	if 'divfaceB' in data: outfile.write("# divergence of face magnetic field (T/m)\n")
	if 'edgeE' in data: outfile.write("# edge-directed electric field on cell edges (V/m, x,-y,-z; x,-y,+z; ...; z,+x,+y)\n")
	if 'faceB' in data: outfile.write("# normal component of magnetic field on faces (T, -x,+x,-y,...,+z)\n")
	if 'mpi' in data: outfile.write("# MPI rank of cell owner\n")
	if 'mhd info' in data: outfile.write("# MHD solver information\n")
	if 'ref lvls' in data: outfile.write("# Minimum and maximum target refinement level\n")
	outfile.write(
		"# Physical constants:\n"
		+ "# adiabatic index: " + str(adiabatic_index)
		+ "\n# vacuum permeability: " + str(vacuum_permeability)
		+ "\n# Simulation time: " + str(sim_time)
		+ "\n# Simulation step: " + str(sim_step) + "\n")

	for i in range(len(cell_ids)):
		cell_id = cell_ids[i]

		# calculate cell geometry, defined by get_center() function in 
		# dccrg_cartesian_geometry.hpp file of dccrg
		cell_id = int(cell_id - 1)
		cell_index = (
			int(cell_id % ref_lvl_0_cells[0]),
			int(cell_id / ref_lvl_0_cells[0] % ref_lvl_0_cells[1]),
			int(cell_id / (ref_lvl_0_cells[0] * ref_lvl_0_cells[1]))
		)

		cell_center = (
			grid_start[0] + lvl_0_cell_length[0] * (0.5 + cell_index[0]),
			grid_start[1] + lvl_0_cell_length[1] * (0.5 + cell_index[1]),
			grid_start[2] + lvl_0_cell_length[2] * (0.5 + cell_index[2])
		)
		cell_id = int(cell_id + 1)

		outfile.write(
			str(cell_center[0]) + ' '
			+ str(cell_center[1]) + ' '
			+ str(cell_center[2]) + ' ')

		if 'mhd' in data:
			infile.seek(data['mhd'] + i*8*8, 0)
			mas, mom, nrj, mag = numpy.fromfile(
				infile,
				dtype = 'double, 3double, double, 3double',
				count = 1
			)[0]

			kin_nrj = (mom[0]**2 + mom[1]**2 + mom[2]**2) / 2 / mas
			mag_nrj = (mag[0]**2 + mag[1]**2 + mag[2]**2) / 2 / vacuum_permeability
			pressure = (nrj - kin_nrj - mag_nrj) * (adiabatic_index - 1)
			outfile.write(
				str(mas/proton_mass) + ' '
				+ str(mom[0]/mas) + ' '
				+ str(mom[1]/mas) + ' '
				+ str(mom[2]/mas) + ' '
				+ str(pressure) + ' '
				+ str(mag[0]) + ' '
				+ str(mag[1]) + ' '
				+ str(mag[2]) + ' ')

		if 'primary' in data:
			infile.seek(data['primary'] + i*(6+12), 0)
			p = numpy.fromfile(infile, dtype = '18u1', count = 1)[0]
			for i in range(len(p)):
				outfile.write(str(p[i]) + ' ')

		if 'bgB' in data:
			infile.seek(data['bgB'] + i*3*6*8, 0)
			B = numpy.fromfile(infile, dtype = '18double', count = 1)[0]
			for i in range(len(B)):
				outfile.write(str(B[i]) + ' ')

		if 'divfaceB' in data:
			infile.seek(data['divfaceB'] + i*8, 0)
			B = numpy.fromfile(infile, dtype = 'double', count = 1)[0]
			outfile.write(str(B) + ' ')

		if 'edgeE' in data:
			infile.seek(data['edgeE'] + i*12*8, 0)
			E = numpy.fromfile(infile, dtype = '12double', count = 1)[0]
			for i in range(len(E)):
				outfile.write(str(E[i]) + ' ')

		if 'faceB' in data:
			infile.seek(data['faceB'] + i*6*8, 0)
			B = numpy.fromfile(infile, dtype = '6double', count = 1)[0]
			for i in range(len(B)):
				outfile.write(str(B[i]) + ' ')

		if 'rank' in data:
			infile.seek(data['rank'] + i*4, 0)
			rank = numpy.fromfile(infile, dtype = 'intc', count = 1)[0]
			outfile.write(str(rank) + ' ')

		if 'mhd info' in data:
			infile.seek(data['mhd info'] + i*4, 0)
			s = numpy.fromfile(infile, dtype = 'uintc', count = 1)[0]
			outfile.write(str(s) + ' ')

		if 'ref lvls' in data:
			infile.seek(data['ref lvls'] + i*8, 0)
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
