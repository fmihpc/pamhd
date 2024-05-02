#! /usr/bin/env python3
'''
Checks differences between MHD outputs of PAMHD.

Copyright 2016 Ilja Honkonen
Copyright 2024 Finnish Meteorogical Institute
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
from os.path import dirname, join
from sys import modules
spec = spec_from_file_location('common', join(dirname(__file__), 'common.py'))
common = module_from_spec(spec)
modules['common'] = common
spec.loader.exec_module(common)


'''
Compares given infiles.

Returns empty string if relative differences are smaller than in args,
otherwise returns description of last difference that was found.
'''
def diff(args, infile1_name, infile2_name, outfile_name):
	from numpy import abs, array, maximum

	try:
		infile1 = open(infile1_name, 'rb')
		infile2 = open(infile2_name, 'rb')
	except Exception as e:
		return "Couldn't open one or more files for reading: " + str(e)

	error = ''

	data1 = common.get_metadata(infile1)
	data2 = common.get_metadata(infile2)

	if data1['file_version'] != data2['file_version']:
		error = 'File version differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['file_version']) + ' != ' + str(data2['file_version'])
	if data1['sim_step'] != data2['sim_step']:
		error = 'Simulation step differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['sim_step']) + ' != ' + str(data2['sim_step'])
	if data1['sim_time'] != data2['sim_time']:
		error = 'Simulation time differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['sim_time']) + ' != ' + str(data2['sim_time'])
	if data1['adiabatic_index'] != data2['adiabatic_index']:
		error = 'Adiabatic index differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['adiabatic_index']) + ' != ' + str(data2['adiabatic_index'])
	if data1['proton_mass'] != data2['proton_mass']:
		error = 'Proton mass differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['proton_mass']) + ' != ' + str(data2['proton_mass'])
	if data1['vacuum_permeability'] != data2['vacuum_permeability']:
		error = 'Vacuum permeability differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['vacuum_permeability']) + ' != ' + str(data2['vacuum_permeability'])
	if data1['endianness'] != data2['endianness']:
		error = 'Endianness differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['endianness']) + ' != ' + str(data2['endianness'])
	if data1['ref_lvl_0_cells'] != data2['ref_lvl_0_cells']:
		error = 'Refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['ref_lvl_0_cells']) + ' != ' + str(data2['ref_lvl_0_cells'])
	if data1['max_ref_lvl'] != data2['max_ref_lvl']:
		error = 'Maximum refinement level differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['max_ref_lvl']) + ' != ' + str(data2['max_ref_lvl'])
	if data1['neighborhood_length'] != data2['neighborhood_length']:
		error = 'Neighborhood length differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['neighborhood_length']) + ' != ' + str(data2['neighborhood_length'])
	if data1['periodicity'] != data2['periodicity']:
		error = 'Grid periodicity differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['periodicity']) + ' != ' + str(data2['periodicity'])
	if data1['geometry_id'] != data2['geometry_id']:
		error = 'Geometry id differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['geometry_id']) + ' != ' + str(data2['geometry_id'])
	if data1['grid_start'] != data2['grid_start']:
		error = 'Grid start coordinate differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['grid_start']) + ' != ' + str(data2['grid_start'])
	if data1['lvl_0_cell_length'] != data2['lvl_0_cell_length']:
		error = 'Length of refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['lvl_0_cell_length']) + ' != ' + str(data2['lvl_0_cell_length'])
	if data1['total_cells'] != data2['total_cells']:
		error = 'Total number of cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['total_cells']) + ' != ' + str(data2['total_cells'])
	data1['cells'] = common.get_cells(infile1, data1['total_cells'])
	data2['cells'] = common.get_cells(infile2, data2['total_cells'])
	ids1 = set([i for i in data1['cells']])
	ids2 = set([i for i in data2['cells']])
	if ids1 != ids2:
		error = 'Cells that exist differs between ' + infile1_name + ' and ' + infile2_name
	if error != '':
		return error

	sim_data1_ = common.get_cell_data(infile1, data1, range(len(data1['cells'])))
	sim_data1 = dict()
	for i in range(len(data1['cells'])):
		sim_data1[data1['cells'][i]] = sim_data1_[i]

	sim_data2_ = common.get_cell_data(infile2, data2, range(len(data2['cells'])))
	sim_data2 = dict()
	for i in range(len(data2['cells'])):
		sim_data2[data2['cells'][i]] = sim_data2_[i]

	for i in range(len(data1['cells'])):
		cell_id = data1['cells'][i]

		rho1 = sim_data1[cell_id]['mhd'][0]
		rho2 = sim_data2[cell_id]['mhd'][0]
		denom = maximum(abs(rho1), abs(rho2))
		if denom != 0 and args.mass < abs(rho1 - rho2) / denom:
			error = 'Mass density differs too much in cell ' + str(cell_id) + ': ' + str(rho1) + ' vs ' + str(rho2)

		vel1 = array(sim_data1[cell_id]['mhd'][1]) / rho1
		vel1 = vel1.dot(vel1)**0.5
		vel2 = array(sim_data2[cell_id]['mhd'][1]) / rho2
		vel2 = vel2.dot(vel2)**0.5
		denom = maximum(vel1, vel2)
		if denom != 0 and args.velocity < abs(vel1 - vel2) / denom:
			error = 'Velocity differs too much in cell ' + str(cell_id) + ': ' + str(vel1) + ' vs ' + str(vel2)

		nrj1 = sim_data1[cell_id]['mhd'][2]
		nrj2 = sim_data2[cell_id]['mhd'][2]
		denom = maximum(abs(nrj1), abs(nrj2))
		if denom != 0 and args.energy < abs(nrj1 - nrj2) / denom:
			error = 'Energy density differs too much in cell ' + str(cell_id) + ': ' + str(nrj1) + ' vs ' + str(nrj2)

		mag1 = array(sim_data1[cell_id]['mhd'][3])
		mag1 = mag1.dot(mag1)**0.5
		mag2 = array(sim_data2[cell_id]['mhd'][3])
		mag2 = mag2.dot(mag2)**0.5
		denom = maximum(mag1, mag2)
		if denom != 0 and args.magnetic_field < abs(mag1 - mag2) / denom:
			error = 'Perturbed magnetic field differs too much in cell ' + str(cell_id) + ': ' + str(mag1) + ' vs ' + str(mag2)

	return error


if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	from os.path import basename
	from sys import argv

	parser = ArgumentParser(
		formatter_class = ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--mass', type = float, default = 1e-9, help = 'Maximum allowed relative difference in mass density')
	parser.add_argument('--velocity', type = float, default = 1e-9, help = 'Maximum allowed relative difference in velocity magnitude')
	parser.add_argument('--energy', type = float, default = 1e-9, help = 'Maximum allowed relative difference in energy density')
	parser.add_argument('--magnetic-field', type = float, default = 1e-9, help = 'Maximum allowed relative difference in perturbed magnetic field magnitude')
	parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of files to compare ([--files] diff1_in1 diff1_in2 diff2_in1 diff2_in2 ...)')
	args = parser.parse_args()

	if len(args.files) % 2 > 0:
		exit('Even number of file arguments required')

	ok = True
	for i in range(0, len(args.files), 2):
		error = diff(args, args.files[i], args.files[i+1],
			join(dirname(args.files[i]), 'diff_' + basename(args.files[i])))
		if error != '':
			ok = False
	if ok:
		exit()
	else:
		exit(error)
