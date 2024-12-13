#! /usr/bin/env python3
'''
Converts MHD output of PAMHD test program to ASCII VTK format.

Copyright 2024 Finnish Meteorological Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


Author(s): Ilja Honkonen
'''

from argparse import ArgumentParser
from imp import load_source
from os.path import basename, dirname, join
from sys import argv, stdout
common = load_source('common', join(dirname(__file__), 'common.py'))


parser = ArgumentParser()
parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of input files to convert to vtk.')
parser.add_argument('--cells', metavar = 'C', type = int, default = -1, help = 'If > 0, process file(s) C cells at a time to decrease required memory.')
args = parser.parse_args()

for filename_ in args.files:
	with open(filename_, 'rb') as infile:
		meta = common.get_metadata(infile)
		if meta['file_version'] != 3:
			print('Unsupported file version:', meta['file_version'])
			continue
		if meta['geometry_id'] != 1:
			print('Unsupported geometry:', meta['geometry_id'])
			continue
		meta['cells'] = common.get_cells(infile, meta['total_cells'])

		with open(filename_.replace('dc', 'vtk'), 'w') as outfile:
			common.write_metadata_vtk(outfile, meta)
			for variable in meta['var_data_start'].keys():
				common.write_variable_vtk(infile, outfile, variable, meta, nr_cells = args.cells)
