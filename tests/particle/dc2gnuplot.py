#! /usr/bin/env python3
'''
Plots output from particle test program of PAMHD using gnuplot.

Copyright 2015, 2016 Ilja Honkonen
Copyright 2024, 2025 Finnish Meteorological Institute
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
from os.path import dirname, isfile, join, realpath
from pathlib import Path
from sys import modules
spec = spec_from_file_location(
	'common',
	join(Path(realpath(__file__)).parent.parent, 'common.py')
)
common = module_from_spec(spec)
modules['common'] = common
spec.loader.exec_module(common)


def plot1d(data, outname):
	rl0c = data['ref_lvl_0_cells']
	dim = None
	if rl0c[0] > 1:
		dim = 0
	elif rl0c[1] > 1:
		dim = 1
	else:
		dim = 2
	common_1d = 'set term png enhanced size 800, 600'
	gs = data['grid_start']
	l0cl = data['lvl_0_cell_length']
	ge = (
		gs[0] + rl0c[0] * l0cl[0],
		gs[1] + rl0c[1] * l0cl[1],
		gs[2] + rl0c[2] * l0cl[2]
	)
	with open(outname + '.gnuplot', 'w') as outfile:
		outfile.write(
			common_1d + '\nset title "Time '
			+ str(data['sim_time']) + '"\nset xlabel "Dimension '
			+ str(dim+1) + '"\nset xrange [' + str(gs[dim])
			+ ':' + str(ge[dim]) + ']\n')

		if 'volE    ' in data['cell_data'][0]:
			outfile.write(
				'set title "Time ' + str(data['sim_time'])
				+ '"\nset xlabel "Dimension ' + str(dim+1)
				+ '"\nset xrange [' + str(gs[dim]) + ':'
				+ str(ge[dim]) + ']\nset output "'
				+ outname + '_vE.png"\n'
				+ 'set ylabel "Volume electric field"\n'
				+ 'set key horizontal bottom outside\nplot '
				+ '"-" u 1:2 lw 2 t "Dimension 1", '
				+ '"-" u 1:2 lw 2 t "Dimension 2", '
				+ '"-" u 1:2 lw 2 t "Dimension 3"\n')
			for i in range(len(data['cells'])):
				outfile.write(
					str(data['centers'][i][dim]) + ' '
					+ str(data['cell_data'][i]['volE    '][0]) + '\n')
			outfile.write('end\n')
			for i in range(len(data['cells'])):
				outfile.write(
					str(data['centers'][i][dim]) + ' '
					+ str(data['cell_data'][i]['volE    '][1]) + '\n')
			outfile.write('end\n')
			for i in range(len(data['cells'])):
				outfile.write(
					str(data['centers'][i][dim]) + ' '
					+ str(data['cell_data'][i]['volE    '][2]) + '\n')
			outfile.write('end\n')


if __name__ == '__main__':
	from argparse import ArgumentParser
	from subprocess import run

	parser = ArgumentParser()
	parser.add_argument('files', metavar = 'F', nargs = '*', help = 'dc files to plot')
	parser.add_argument('--verbose', action = 'store_true', default = False)
	args = parser.parse_args()

	if args.verbose:
		stdout.write('Plotting files: ' + str(args.files) + '\n')

	for inname in args.files:
		with open(inname, 'rb') as infile:
			data = common.get_metadata(infile)
			data['cells'] = common.get_cells(infile, data['total_cells'])
			if len(data['cells']) == 0:
				continue
			data['cell_data'] = common.get_cell_data(infile, data, range(len(data['cells'])))

			if data['file_version'] != 4:
				exit('Unsupported file version: ' + str(data['file_version']))

			if data['geometry_id'] != 1:
				exit('Unsupported geometry: ' + str(data['geometry_id']))

			data['centers'] = []
			for i in range(len(data['cells'])):
				cell_id = data['cells'][i]
				center, _ = common.get_cell_geom(data, cell_id)
				data['centers'].append(center)

			rl0c = data['ref_lvl_0_cells']
			dims = 0
			for c in rl0c:
				if c > 1:
					dims += 1
			if dims == 0:
				continue
			if dims > 1:
				print('Only 1d plot supported')
				continue

			plot1d(data, inname[:-3])
			run('gnuplot ' + inname[:-3] + '.gnuplot', shell = True)
