#! /usr/bin/env python3
'''
Extracts subvolumes of MHD output of PAMHD.

Copyright 2016 Ilja Honkonen
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


from imp import load_source
from os.path import basename, dirname, join
common = load_source('common', join(dirname(__file__), 'common.py'))

'''
Reads simulation data from infile_name and writes it into outfile_name while excluding cells whose center is outside of given volume
'''
def extract(
	infile_name,
	outfile_name,
	min_x = float('-inf'),
	max_x = float('+inf'),
	min_y = float('-inf'),
	max_y = float('+inf'),
	min_z = float('-inf'),
	max_z = float('+inf')
):
	from os.path import exists, isfile
	from numpy import dtype, fromfile, int32, uint64, zeros

	try:
		infile = open(infile_name, 'rb')
	except Exception as e:
		print("Couldn't open/read from", infile_name, ':', e)
		return

	try:
		meta = common.get_metadata(infile)
	except Exception as e:
		print('Failed to load metadata from', infile_name, ':', e)
		return

	if meta['file_version'] != 3:
		print('File version', meta['file_version'], 'not supported:', infile_name)
		return

	try:
		incells = common.get_cells(infile, meta['total_cells'])
	except Exception as e:
		print("Couldn't load cell list from", infile_name, ':', e)
		return

	# only read&write cells overlapping with requested volume
	outcells = []
	inds = [] # location of outcell in incells
	for ind in range(len(incells)):
		cell = incells[ind]
		center, length = common.get_cell_geom(meta, cell)
		cix = center[0] - length[0]/2 # cell_min_x
		cax = center[0] + length[0]/2 # cell_max_x
		ciy = center[1] - length[1]/2
		cay = center[1] + length[1]/2
		ciz = center[2] - length[2]/2
		caz = center[2] + length[2]/2
		if cax >= min_x and cix <= max_x \
			and cay >= min_y and ciy <= max_y \
			and caz >= min_z and ciz <= max_z \
		:
			outcells.append(cell)
			inds.append(ind)

	try:
		data = common.get_cell_data(infile, meta, inds)
	except Exception as e:
		print("Couldn't load simulation data from", infile_name, ':', e)
		return

	try:
		outfile = open(outfile_name, 'wb')
	except Exception as e:
		print("Couldn't open/write to", outfile_name, ':', e)
		return

	common.set_metadata(outfile, meta)
	common.set_cells(outfile, outcells)
	saved_vars = common.set_cell_data(outfile, meta, data)
	# write header again with correct values
	meta['total_cells'] = len(outcells)
	nr_vars_orig = len(meta['variable_offsets'])
	meta['variable_offsets'] = [offset[1] for offset in saved_vars]
	while len(meta['variable_offsets']) < nr_vars_orig:
		meta['variable_offsets'].append(meta['variable_offsets'][-1])
	outfile.seek(0, 0)
	common.set_metadata(outfile, meta)


if __name__ == '__main__':
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	from sys import argv

	parser = ArgumentParser(
		formatter_class = ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of input and output files to process (in1 out1 in2 out2 ...)')
	parser.add_argument('--min-x', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center x coordinate larger than MIN_X')
	parser.add_argument('--max-x', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center x coordinate smaller than MAX_X')
	parser.add_argument('--min-y', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center y coordinate larger than MIN_Y')
	parser.add_argument('--max-y', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center y coordinate smaller than MAX_Y')
	parser.add_argument('--min-z', type = float, default = float('-inf'), help = 'Include in output file(s) only cells with center z coordinate larger than MIN_Z')
	parser.add_argument('--max-z', type = float, default = float('+inf'), help = 'Include in output file(s) only cells with center z coordinate smaller than MAX_Z')
	args = parser.parse_args()

	if len(args.files) % 2 > 0:
		exit('Even number of file arguments required (in1 out1 in2 out2 ...)')

	if args.min_x >= args.max_x:
		exit('Minimum x coordinate must be less than maximum')
	if args.min_y >= args.max_y:
		exit('Minimum y coordinate must be less than maximum')
	if args.min_z >= args.max_z:
		exit('Minimum z coordinate must be less than maximum')

	for i in range(0, len(args.files), 2):
		extract(
			args.files[i],
			args.files[i + 1],
			min_x = args.min_x,
			max_x = args.max_x,
			min_y = args.min_y,
			max_y = args.max_y,
			min_z = args.min_z,
			max_z = args.max_z
		)
