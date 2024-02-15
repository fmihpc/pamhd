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


'''
Returns dict of simulation data in given file open for reading
in binary mode, excluding user data stored in simulation cells.
'''
def get_sim_data(infile):
	from numpy import dtype, fromfile, int32, uint64, zeros
	ret_val = dict()
	file_version = fromfile(infile, dtype = 'uint64', count = 1)
	if file_version[0] != 2:
		exit('Unsupported file version: ' + str(file_version[0]))
	ret_val['file_version'] = file_version
	ret_val['sim_step'] = fromfile(infile, dtype = 'uint64', count = 1)
	ret_val['sim_params'] = fromfile(infile, dtype = '4double', count = 1)
	endianness = fromfile(infile, dtype = 'uint64', count = 1)
	if endianness[0] != 0x1234567890abcdef:
		exit('Unsupported endianness: ' + str(endianness[0]))
	ret_val['endianness'] = endianness
	ret_val['ref_lvl_0_cells'] = fromfile(infile, dtype = '3uint64', count = 1)
	ret_val['max_ref_lvl'] = fromfile(infile, dtype = 'intc', count = 1)
	ret_val['neighborhood_length'] = fromfile(infile, dtype = 'uintc', count = 1)
	ret_val['periodicity'] = fromfile(infile, dtype = '3uint8', count = 1)
	geometry_id = fromfile(infile, dtype = 'intc', count = 1)
	if geometry_id[0] != int32(1):
		exit('Unsupported geometry: ' + str(geometry_id[0]))
	ret_val['geometry_id'] = geometry_id
	ret_val['grid_start'] = fromfile(infile, dtype = '3double', count = 1)
	ret_val['lvl_0_cell_length'] = fromfile(infile, dtype = '3double', count = 1)
	total_cells = fromfile(infile, dtype = 'uint64', count = 1)
	ret_val['total_cells'] = total_cells
	ret_val['ids_offsets'] = fromfile(infile, dtype = '2uint64', count = total_cells[0])
	return ret_val


'''
Compares given infiles.

Returns empty string if relative differences are smaller than in args,
otherwise returns description of last difference that was found.

Optionally writes difference between simulation data in
infile1_name and infile2_name into outfile_name.
'''
def diff(args, infile1_name, infile2_name, outfile_name):
	from os.path import exists, isfile
	from numpy import abs, dtype, fromfile, int32, maximum, uint64, zeros

	if not exists(infile1_name):
		return 'Path ' + infile1_name + " doesn't exist"
	if not isfile(infile1_name):
		return 'Path ' + infile1_name + ' is not a file'
	if not exists(infile2_name):
		return 'Path ' + infile2_name + " doesn't exist"
	if not isfile(infile2_name):
		return 'Path ' + infile2_name + ' is not a file'

	try:
		infile1 = open(infile1_name, 'rb')
		infile2 = open(infile2_name, 'rb')
	except Exception as e:
		return "Couldn't open one or more files for reading: " + str(e)

	error = ''

	data1 = get_sim_data(infile1)
	data2 = get_sim_data(infile2)

	if data1['file_version'] != data2['file_version']:
		error = 'File version differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['file_version']) + ' != ' + str(data2['file_version'])
	if data1['sim_step'] != data2['sim_step']:
		error = 'Simulations step differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['sim_step']) + ' != ' + str(data2['sim_step'])
	if (data1['sim_params'] != data2['sim_params']).any():
		error = 'Simulation parameters differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['sim_params']) + ' != ' + str(data2['sim_params'])
	if data1['endianness'] != data2['endianness']:
		error = 'Endianness differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['endianness']) + ' != ' + str(data2['endianness'])
	if (data1['ref_lvl_0_cells'] != data2['ref_lvl_0_cells']).any():
		error = 'Refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['ref_lvl_0_cells']) + ' != ' + str(data2['ref_lvl_0_cells'])
	if data1['max_ref_lvl'] != data2['max_ref_lvl']:
		error = 'Maximum refinement level differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['max_ref_lvl']) + ' != ' + str(data2['max_ref_lvl'])
	if (data1['neighborhood_length'] != data2['neighborhood_length']).any():
		error = 'Neighborhood length differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['neighborhood_length']) + ' != ' + str(data2['neighborhood_length'])
	if (data1['periodicity'] != data2['periodicity']).any():
		error = 'Grid periodicity differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['periodicity']) + ' != ' + str(data2['periodicity'])
	if data1['geometry_id'] != data2['geometry_id']:
		error = 'Geometry id differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['geometry_id']) + ' != ' + str(data2['geometry_id'])
	if (data1['grid_start'] != data2['grid_start']).any():
		error = 'Grid start coordinate differs between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['grid_start']) + ' != ' + str(data2['grid_start'])
	if (data1['lvl_0_cell_length'] != data2['lvl_0_cell_length']).any():
		error = 'Length of refinement level 0 cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['lvl_0_cell_length']) + ' != ' + str(data2['lvl_0_cell_length'])
	if data1['total_cells'] != data2['total_cells']:
		error = 'Total number of cells differ between ' + infile1_name + ' and ' + infile2_name + ': ' + str(data1['total_cells']) + ' != ' + str(data2['total_cells'])
	ids1 = set([i[0] for i in data1['ids_offsets']])
	ids2 = set([i[0] for i in data2['ids_offsets']])
	if ids1 != ids2:
		error = 'Cells that exist differs between ' + infile1_name + ' and ' + infile2_name

	if args.write:
		try:
			outfile = open(outfile_name, 'wb')
		except Exception as e:
			return "Couldn't open output file " + output_file_name + ': ' + str(e)
		data1['file_version'].tofile(outfile)
		data1['sim_step'].tofile(outfile)
		data1['sim_params'].tofile(outfile)
		data1['endianness'].tofile(outfile)
		data1['ref_lvl_0_cells'].tofile(outfile)
		data1['max_ref_lvl'].tofile(outfile)
		data1['neighborhood_length'].tofile(outfile)
		data1['periodicity'].tofile(outfile)
		data1['geometry_id'].tofile(outfile)
		data1['grid_start'].tofile(outfile)
		data1['lvl_0_cell_length'].tofile(outfile)
		data1['total_cells'].tofile(outfile)
		data1['ids_offsets'].tofile(outfile)

	cell_data_t = 'double, 3double, double, 3double, intc, intc, 3double, 3double, 12double, 3double, 3double, 3double, 3double, 3double, 3double, double, intc'
	ids_offsets2 = dict() # cells in different order than in infile1
	for i in data2['ids_offsets']:
		ids_offsets2[i[0]] = i[1]
	for i in range(len(data1['ids_offsets'])):
		cell_id = data1['ids_offsets'][i][0]
		offset = data1['ids_offsets'][i][1]
		infile1.seek(offset, 0)
		cell_data1 = fromfile(infile1, dtype = cell_data_t, count = 1)
		infile2.seek(ids_offsets2[cell_id], 0)
		cell_data2 = fromfile(infile2, dtype = cell_data_t, count = 1)

		rho1 = cell_data1[0][0]
		rho2 = cell_data2[0][0]
		denom = maximum(abs(rho1), abs(rho2))
		if denom != 0 and args.mass < abs(rho1 - rho2) / denom:
			error = 'Mass density differs too much in cell ' + str(cell_id) + ': ' + str(rho1) + ' vs ' + str(rho2)

		vel1 = cell_data1[0][1] / cell_data1[0][0]
		vel1 = vel1.dot(vel1)**0.5
		vel2 = cell_data2[0][1] / cell_data2[0][0]
		vel2 = vel2.dot(vel2)**0.5
		if denom != 0 and args.velocity < abs(vel1 - vel2) / maximum(vel1, vel2):
			error = 'Velocity differs too much in cell ' + str(cell_id) + ': ' + str(vel1) + ' vs ' + str(vel2)

		nrj1 = cell_data1[0][2]
		nrj2 = cell_data2[0][2]
		denom = maximum(abs(nrj1), abs(nrj2))
		if denom != 0 and args.energy < abs(nrj1 - nrj2) / denom:
			error = 'Energy density differs too much in cell ' + str(cell_id) + ': ' + str(nrj1) + ' vs ' + str(nrj2)

		mag1 = cell_data1[0][3]
		mag1 = mag1.dot(mag1)**0.5
		mag2 = cell_data2[0][3]
		mag2 = mag2.dot(mag2)**0.5
		if denom != 0 and args.magnetic_field < abs(mag1 - mag2) / maximum(mag1, mag2):
			error = 'Perturbed magnetic field differs too much in cell ' + str(cell_id) + ': ' + str(mag1) + ' vs ' + str(mag2)

		if args.write:
			out_cell_data = cell_data1.copy()
			for i in range(len(out_cell_data[0])):
				out_cell_data[0][i] = abs(cell_data1[0][i] - cell_data2[0][i])
			out_cell_data.tofile(outfile)

	return error


if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	from os.path import basename, dirname, join
	from sys import argv

	parser = ArgumentParser(
		formatter_class = ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--write', action = 'store_true', help = 'Write differences between every 2 input files into output file with name based on first input file')
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
