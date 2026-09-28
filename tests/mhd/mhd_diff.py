#! /usr/bin/env python3
'''
Checks differences between MHD outputs of PAMHD.

Copyright 2016 Ilja Honkonen
Copyright 2024, 2025, 2026 Finnish Meteorogical Institute
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

from copy import deepcopy
from importlib.util import module_from_spec, spec_from_file_location
from os.path import dirname, join, realpath
from pathlib import Path
from sys import modules
spec = spec_from_file_location(
	'common',
	join(Path(realpath(__file__)).parent.parent, 'common.py')
)
common = module_from_spec(spec)
modules['common'] = common
spec.loader.exec_module(common)


def transform_meta(args, data):
	if args.rotate == None:
		return data

	i1, i2 = -1, -1
	if args.rotate == 'x':
		i1, i2 = 1, 2
	if args.rotate == 'y':
		i1, i2 = 0, 2
	if args.rotate == 'z':
		i1, i2 = 0, 1

	temp = data['ref_lvl_0_cells'][i1]
	data['ref_lvl_0_cells'][i1] = data['ref_lvl_0_cells'][i2]
	data['ref_lvl_0_cells'][i2] = temp
	temp = data['periodicity'][i1]
	data['periodicity'][i1] = data['periodicity'][i2]
	data['periodicity'][i2] = temp
	temp = data['grid_start'][i1]
	data['grid_start'][i1] = data['grid_start'][i2]
	data['grid_start'][i2] = temp
	temp = data['lvl_0_cell_length'][i1]
	data['lvl_0_cell_length'][i1] = data['lvl_0_cell_length'][i2]
	data['lvl_0_cell_length'][i2] = temp
	return data


def transform_cells(args, old_meta, new_meta, cells):
	if args.rotate == None:
		return cells

	rl0c = old_meta['ref_lvl_0_cells']
	mrl = old_meta['max_ref_lvl']
	max_index = (
		rl0c[0]*2**mrl - 1,
		rl0c[1]*2**mrl - 1,
		rl0c[2]*2**mrl - 1)
	new_cells = []
	for old_cell in cells:
		rlvl, old_index = common.get_cell_info(old_meta, old_cell)
		if args.rotate == 'x':
			new_index = (
				old_index[0],
				max_index[2] - old_index[2],
				old_index[1])
		elif args.rotate == 'y':
			new_index = (
				max_index[2] - old_index[2],
				old_index[1],
				old_index[0])
		elif args.rotate == 'z':
			new_index = (
				max_index[1] - old_index[1],
				old_index[0],
				old_index[2])
		new_cell = common.get_cell(new_meta, rlvl, new_index)
		new_cells.append(new_cell)
	return new_cells


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
		raise RuntimeError("Couldn't open " + infile1_name + ' or ' + infile2_name + ' for reading: ' + repr(e))

	data1 = common.get_metadata(infile1)
	orig2 = common.get_metadata(infile2)
	data2 = transform_meta(args, deepcopy(orig2))

	if data1['file_version'] != data2['file_version']:
		raise RuntimeError('File version differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['file_version']) + ' != ' + str(data2['file_version']))

	sim_step1 = data1['sim_step']
	sim_step2 = data2['sim_step']
	denom = maximum(sim_step1, sim_step2)
	diff_ = abs(sim_step1 - sim_step2) - args.step_min
	if denom > 0 and args.step < diff_ / denom:
		raise RuntimeError('Simulation step differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(sim_step1) + ' != ' + str(sim_step2))

	sim_time1 = data1['sim_time']
	sim_time2 = data2['sim_time']
	denom = maximum(abs(sim_time1), abs(sim_time2))
	diff_ = abs(sim_time1 - sim_time2) - args.time_min
	if denom > 0 and args.time < diff_ / denom:
		raise RuntimeError('Simulation time differs too much between ' + infile1_name + ' and ' + infile2_name + ': ' + str(abs(sim_time1 - sim_time2)) + ' (' + str(sim_time1) + ' vs ' + str(sim_time2) + ')')

	if data1['adiabatic_index'] != data2['adiabatic_index']:
		raise RuntimeError('Adiabatic index differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['adiabatic_index']) + ' != ' + str(data2['adiabatic_index']))
	if data1['proton_mass'] != data2['proton_mass']:
		raise RuntimeError('Proton mass differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['proton_mass']) + ' != ' + str(data2['proton_mass']))
	if data1['vacuum_permeability'] != data2['vacuum_permeability']:
		raise RuntimeError('Vacuum permeability differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['vacuum_permeability']) + ' != ' + str(data2['vacuum_permeability']))
	if data1['endianness'] != data2['endianness']:
		raise RuntimeError('Endianness differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['endianness']) + ' != ' + str(data2['endianness']))

	rl0c1 = data1['ref_lvl_0_cells']
	rl0c2 = data2['ref_lvl_0_cells']
	if not args.cells and rl0c1 != rl0c2:
		raise RuntimeError('Refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(rl0c1) + ' != ' + str(rl0c2))
	if not args.cells and data1['max_ref_lvl'] != data2['max_ref_lvl']:
		raise RuntimeError('Maximum refinement level differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['max_ref_lvl']) + ' != ' + str(data2['max_ref_lvl']))
	if data1['neighborhood_length'] != data2['neighborhood_length']:
		raise RuntimeError('Neighborhood length differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['neighborhood_length']) + ' != ' + str(data2['neighborhood_length']))
	if data1['periodicity'] != data2['periodicity']:
		raise RuntimeError('Grid periodicity differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['periodicity']) + ' != ' + str(data2['periodicity']))
	if data1['geometry_id'] != data2['geometry_id']:
		raise RuntimeError('Geometry id differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['geometry_id']) + ' != ' + str(data2['geometry_id']))
	if data1['grid_start'] != data2['grid_start']:
		raise RuntimeError('Grid start coordinate differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['grid_start']) + ' != ' + str(data2['grid_start']))
	if not args.cells and data1['lvl_0_cell_length'] != data2['lvl_0_cell_length']:
		raise RuntimeError('Length of refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['lvl_0_cell_length']) + ' != ' + str(data2['lvl_0_cell_length']))
	if not args.cells and data1['total_cells'] != data2['total_cells']:
		raise RuntimeError('Total number of cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['total_cells']) + ' != ' + str(data2['total_cells']))
	cells1 = common.get_cells(infile1, data1['total_cells'])
	orig_cells2 = common.get_cells(infile2, orig2['total_cells'])
	cells2 = transform_cells(args, orig2, data2, orig_cells2)

	ids1 = set([c for c in cells1])
	ids2 = set([c for c in cells2])
	if not args.cells and ids1 != ids2:
		raise RuntimeError('Cells that exist differs between ' + infile1_name + ' and ' + infile2_name)

	cell2orig = dict()
	for i in range(len(cells2)):
		cell2orig[cells2[i]] = orig_cells2[i]

	# assuming *_in1 files use better resolution, 1 index of *_in2 spans this many:
	index_factor = [1, 1, 1]
	if rl0c1 != rl0c2:
		index_factor = (
			rl0c1[0] // rl0c2[0],
			rl0c1[1] // rl0c2[1],
			rl0c1[2] // rl0c2[2])
		if index_factor[0] <= 0 or index_factor[1] <= 0 or index_factor[2] <= 0:
			print(index_factor[0], index_factor[1], index_factor[2])
			raise RuntimeError('Invalid index scaling factor(s) between ' + infile1_name + ' and ' + infile2_name)

	max_rho, max_nrj, max_vel, max_mag = 0, 0, 0, 0
	for cell1 in cells1:
		if args.cells:
			# corresponding larger cell
			ref_lvl1, index1 = common.get_cell_info(data1, cell1)
			index2 = (
				index1[0] // index_factor[0],
				index1[1] // index_factor[1],
				index1[2] // index_factor[2])
			for ref_lvl2 in range(data2['max_ref_lvl'] + 1):
				cell2 = common.get_cell(data2, ref_lvl2, index2)
				if cell2 in ids2:
					break
		else:
			cell2 = cell2orig[cell1]

		try:
			cell1_pos = cells1.index(cell1)
			cell2_pos = orig_cells2.index(cell2)
		except Exception as e:
			print("Couldn't locate cell(s) in cell list(s) of files", infile1_name, infile2_name)
			raise RuntimeError('qwerty')
		sim_data1 = common.get_cell_data(infile1, data1, range(cell1_pos, cell1_pos + 1), ['mhd     '])[0]['mhd     ']
		sim_data2 = common.get_cell_data(infile2, orig2, range(cell2_pos, cell2_pos + 1), ['mhd     '])[0]['mhd     ']

		rho1 = sim_data1[0]
		rho2 = sim_data2[0]
		if args.norm:
			max_rho = maximum(max_rho, abs(rho1 - rho2))
		else:
			denom = maximum(abs(rho1), abs(rho2))
			diff_ = abs(rho1 - rho2) - args.mass_min
			if denom != 0 and args.mass < diff_ / denom:
				raise RuntimeError('Mass density differs too much in cell ' + str(cell1) + ' of file ' + infile1_name + ' and cell ' + str(cell2) + ' of file ' + infile2_name + ': ' + str(abs(rho1 - rho2)) + ' (' + str(rho1) + ' vs ' + str(rho2) + ')')

		vel1 = array(sim_data1[1]) / rho1
		vel1 = vel1.dot(vel1)**0.5
		vel2 = array(sim_data2[1]) / rho2
		vel2 = vel2.dot(vel2)**0.5
		if args.norm:
			max_vel = maximum(max_vel, abs(vel1 - vel2))
		else:
			denom = maximum(vel1, vel2)
			diff_ = abs(vel1 - vel2) - args.velocity_min
			if denom != 0 and args.velocity < diff_ / denom:
				raise RuntimeError('Velocity differs too much in cell ' + str(cell1) + ' of file ' + infile1_name + ' and cell ' + str(cell2) + ' of file ' + infile2_name + ': ' + str(abs(vel1 - vel2)) + ' (' + str(vel1) + ' vs ' + str(vel2) + ')')

		nrj1 = sim_data1[2]
		nrj2 = sim_data2[2]
		if args.norm:
			max_nrj = maximum(max_nrj, abs(nrj1 - nrj2))
		else:
			denom = maximum(abs(nrj1), abs(nrj2))
			diff_ = abs(nrj1 - nrj2) - args.energy_min
			if denom != 0 and args.energy < diff_ / denom:
				raise RuntimeError('Energy density differs too much in cell ' + str(cell1) + ' of file ' + infile1_name + ' and cell ' + str(cell2) + ' of file ' + infile2_name + ': ' + str(abs(nrj1 - nrj2)) + ' (' + str(nrj1) + ' vs ' + str(nrj2) + ')')

		mag1 = array(sim_data1[3])
		mag1 = mag1.dot(mag1)**0.5
		mag2 = array(sim_data2[3])
		mag2 = mag2.dot(mag2)**0.5
		if args.norm:
			max_mag = maximum(max_mag, abs(mag1 - mag2))
		else:
			denom = maximum(mag1, mag2)
			diff_ = abs(mag1 - mag2) - args.magnetic_field_min
			if denom != 0 and args.magnetic_field < diff_ / denom:
				raise RuntimeError('Perturbed magnetic field differs too much in cell ' + str(cell1) + ' of file ' + infile1_name + ' and cell ' + str(cell2) + ' of file ' + infile2_name + ': ' + str(abs(mag1 - mag2)) + ' (' + str(mag1) + ' vs ' + str(mag2) + ')')

	if args.norm:
		print('#', infile1_name, infile2_name)
		print(max_rho, max_nrj, max_vel, max_mag, flush = True)


if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	from os.path import basename
	from sys import argv

	parser = ArgumentParser(
		formatter_class = ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--cells', default = False, action = 'store_true',
		help = 'Allow differences in number of cells and their logical sizes (use larger file(s) as *_in1)')
	parser.add_argument('--norm', default = False, action = 'store_true',
		help = 'Print maximum norm between MHD variables instead of checking for differences')
	parser.add_argument('--time', type = float, default = 1e-15,
		help = 'Maximum allowed relative difference in simulation time')
	parser.add_argument('--time-min', type = float, default = 1e-15,
		help = 'Minimum value of difference in simulation time')
	parser.add_argument('--step', type = int, default = 0,
		help = 'Maximum allowed relative difference in simulation step')
	parser.add_argument('--step-min', type = int, default = 0,
		help = 'Minimum value of difference in simulation step')
	parser.add_argument('--mass', type = float, default = 1e-15,
		help = 'Maximum allowed relative difference in mass density')
	parser.add_argument('--mass-min', type = float, default = 1e-15,
		help = 'Minimum value of difference in mass density')
	parser.add_argument('--velocity', type = float, default = 1e-15,
		help = 'Maximum allowed relative difference in velocity magnitude')
	parser.add_argument('--velocity-min', type = float, default = 1e-15,
		help = 'Value substracted from velocity difference')
	parser.add_argument('--energy', type = float, default = 1e-15,
		help = 'Maximum allowed relative difference in energy density')
	parser.add_argument('--energy-min', type = float, default = 1e-15,
		help = 'Minimum value of difference in energy density')
	parser.add_argument('--magnetic-field', type = float, default = 1e-15,
		help = 'Maximum allowed relative difference in perturbed magnetic field magnitude')
	parser.add_argument('--magnetic-field-min', type = float, default = 1e-15,
		help = 'Minimum value of difference in magnetic field')
	parser.add_argument('--rotate', metavar = 'A', type = str, default = None,
		help = 'File(s) *_in2 are rotated around axis A before comparing to *_in1')
	parser.add_argument('files', metavar = 'F', nargs = '*',
		help = 'Names of files to compare ([--files] diff1_in1 diff1_in2 diff2_in1 diff2_in2 ...)')
	args = parser.parse_args()

	if len(args.files) % 2 > 0:
		exit('Even number of file arguments required')

	if args.time < 0:
		exit('--time cannot be < 0')
	if args.time_min < 0:
		exit('--time-min cannot be < 0')
	if args.step < 0:
		exit('--step cannot be < 0')
	if args.step_min < 0:
		exit('--step-min cannot be < 0')
	if args.mass < 0:
		exit('--mass cannot be < 0')
	if args.mass_min < 0:
		exit('--mass-min cannot be < 0')
	if args.velocity < 0:
		exit('--velocity cannot be < 0')
	if args.velocity_min < 0:
		exit('--velocity-min cannot be < 0')
	if args.energy < 0:
		exit('--energy cannot be < 0')
	if args.energy_min < 0:
		exit('--energy-min cannot be < 0')
	if args.magnetic_field < 0:
		exit('--magnetic-field cannot be < 0')
	if args.magnetic_field_min < 0:
		exit('--magnetic-field-min cannot be < 0')
	if not args.rotate in {None, 'x', 'y', 'z'}:
		exit('If given, --rotate must be x, y or z')

	if args.norm:
		print('# Maximum norm in density, energy, velocity and magnetic field magnitude')

	ok = True
	for i in range(0, len(args.files), 2):
		try:
			diff(args, args.files[i], args.files[i+1],
				join(dirname(args.files[i]),
				'diff_' + basename(args.files[i])))
		except Exception as e:
			ok = False
			print(repr(e))
	if ok:
		exit()
	else:
		exit(1)
