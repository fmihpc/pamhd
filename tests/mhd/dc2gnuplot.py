#! /usr/bin/env python3
'''
Plots output from MHD test program of PAMHD using gnuplot.

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
from math import isnan
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
		cells = data['cells']
		# skip dont_solve cells
		if len(cells) > 0 and 'mhd info' in data[cells[0]]:
			cells = [cell for cell in data['cells'] if data[cell]['mhd info'] >= 0]

		outfile.write(
			common_1d + '\nset title "Time '
			+ str(data['sim_time']) + '"\nset xlabel "Dimension '
			+ str(dim+1) + '"\nset xrange [' + str(gs[dim])
			+ ':' + str(ge[dim]) + ']\n')

		if 'mhd     ' in data[cells[0]]:
			outfile.write(
				'set output "' + outname + '_n.png"\n'
				+ 'set ylabel "Number density" textcolor lt 1\n'
				+ 'set y2label "Pressure" textcolor lt 3\n'
				+ 'unset key\nset ytics nomirror\nset format x "%.2e"\n'
				+ 'set format y "%.2e"\nset format y2 "%.2e"\n'
				+ 'set y2tics auto\nplot "-" using 1:2 axes '
				+  'x1y1 linewidth 2, "-" u 1:2 axes x1y2 lt 3 lw 2\n')
			for cell in cells:
				mas = data[cell]['mhd     '][0]
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(mas / data['proton_mass']) + '\n')
			outfile.write('end\n')
			for cell in cells:
				mas = data[cell]['mhd     '][0]
				mom = data[cell]['mhd     '][1]
				nrj = data[cell]['mhd     '][2]
				mag = data[cell]['mhd     '][3]
				kin_nrj = (mom[0]**2 + mom[1]**2 + mom[2]**2) / 2 / mas
				if isnan(kin_nrj):
					kin_nrj = 0
				mag_nrj = (mag[0]**2 + mag[1]**2 + mag[2]**2) / 2 / data['vacuum_permeability']
				pressure = (nrj - kin_nrj - mag_nrj) * (data['adiabatic_index'] - 1)
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(pressure) + '\n')
			outfile.write('end\nreset\n')
			outfile.write(
				'set title "Time ' + str(data['sim_time'])
				+ '"\nset xlabel "Dimension ' + str(dim+1)
				+ '"\nset xrange [' + str(gs[dim]) + ':'
				+ str(ge[dim]) + ']\nset output "'
				+ outname + '_V.png"\nset ylabel "Velocity"\n'
				+ 'set key horizontal bottom outside\nplot '
				+ '"-" u 1:2 lw 2 t "component 1", '
				+ '"-" u 1:2 lw 2 t "component 2", '
				+ '"-" u 1:2 lw 2 t "component 3"\n')
			for cell in cells:
				mas = data[cell]['mhd     '][0]
				mom = data[cell]['mhd     '][1]
				vx = mom[0] / mas
				if isnan(vx):
					vx = 0
				outfile.write(str(data['centers'][cell][dim])+' '+str(vx)+'\n')
			outfile.write('end\n')
			for cell in cells:
				mas = data[cell]['mhd     '][0]
				mom = data[cell]['mhd     '][1]
				vy = mom[1] / mas
				if isnan(vy):
					vy = 0
				outfile.write(str(data['centers'][cell][dim])+' '+str(vy)+'\n')
			outfile.write('end\n')
			for cell in cells:
				mas = data[cell]['mhd     '][0]
				mom = data[cell]['mhd     '][1]
				vz = mom[2] / mas
				if isnan(vz):
					vz = 0
				outfile.write(str(data['centers'][cell][dim])+' '+str(vz)+'\n')
			outfile.write('end\nreset\n')
			outfile.write(
				'set title "Time ' + str(data['sim_time'])
				+ '"\nset xlabel "Dimension ' + str(dim+1)
				+ '"\nset xrange [' + str(gs[dim]) + ':'
				+ str(ge[dim]) + ']\nset output "'
				+ outname + '_vB1.png"\nset ylabel "Volume B1"\n'
				+ 'set key horizontal bottom outside\nplot '
				+ '"-" u 1:2 lw 2 t "component 1", '
				+ '"-" u 1:2 lw 2 t "component 2", '
				+ '"-" u 1:2 lw 2 t "component 3"\n')
			for cell in cells:
				mag = data[cell]['mhd     '][3]
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(mag[0]) + '\n')
			outfile.write('end\n')
			for cell in cells:
				mag = data[cell]['mhd     '][3]
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(mag[1]) + '\n')
			outfile.write('end\n')
			for cell in cells:
				mag = data[cell]['mhd     '][3]
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(mag[2]) + '\n')
			outfile.write('end\nreset\n')

		if 'substmin' in data[cells[0]]:
			outfile.write(
				'set title "Time ' + str(data['sim_time'])
				+ '"\nset xlabel "Dimension ' + str(dim+1)
				+ '"\nset xrange [' + str(gs[dim]) + ':'
				+ str(ge[dim]) + ']\nset output "'
				+ outname + '_substep.png"\n'
				+ 'set ylabel "Substepping period"\n')
			if 'timestep' in data[cells[0]]:
				outfile.write(
					'set y2label "Cell dt"\n'
					+ 'set y2tics\nset ytics nomirror\n')
			outfile.write(
				'set key horizontal bottom outside\nplot '
				+ '"-" u 1:2 lw 2 t "min substep", '
				+ '"-" u 1:2 lw 2 t "max substep", '
				+ '"-" u 1:2 lw 2 t "substep"')
			if 'timestep' in data[cells[0]]:
				outfile.write(', "-" u 1:2 axis x1y2 lw 2 t "dt"\n')
			else:
				outfile.write('\n')
			for i in range(len(data['cells'])):
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					# min and max are kept in 2^N format
					+ str(2**data[cell]['substmin']) + '\n')
			outfile.write('end\n')
			for cell in cells:
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(2**data[cell]['substmax']) + '\n')
			outfile.write('end\n')
			for cell in cells:
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(data[cell]['substep ']) + '\n')
			outfile.write('end\n')
		if 'timestep' in data[cells[0]]:
			for cell in cells:
				outfile.write(
					str(data['centers'][cell][dim]) + ' '
					+ str(data[cell]['timestep']) + '\n')
			outfile.write('end\nreset\n')

		if 'bgB     ' in data[cells[0]]:
			outfile.write(
				'set title "Time ' + str(data['sim_time'])
				+ '"\nset xlabel "Dimension ' + str(dim+1)
				+ '"\nset xrange [' + str(gs[dim]) + ':'
				+ str(ge[dim]) + ']\nset output "'
				+ outname + '_B0.png"\nset ylabel "B0"\n'
				+ 'set key horizontal bottom outside\nplot '
				+ '"-" u 1:2 lw 2 t "component 1", '
				+ '"-" u 1:2 lw 2 t "component 2", '
				+ '"-" u 1:2 lw 2 t "component 3"\n')
			for bgB_i in [3, 10, 17]:
				for cell in cells:
					bgB = data[cell]['bgB     ']
					outfile.write(
						str(data['centers'][cell][dim]) + ' '
						+ str(bgB[bgB_i]) + '\n')
				outfile.write('end\n')
			outfile.write('reset\n')


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
			if data['file_version'] != 4:
				exit('Unsupported file version: ' + str(data['file_version']))

			if data['geometry_id'] != 1:
				exit('Unsupported geometry: ' + str(data['geometry_id']))

			data['cells'] = common.get_cells(infile, data['total_cells'])
			if len(data['cells']) == 0:
				continue
			sim_data_ = common.get_cell_data(infile, data, range(len(data['cells'])))
			for i in range(len(data['cells'])):
				data[data['cells'][i]] = sim_data_[i]

			data['centers'] = dict()
			for i in range(len(data['cells'])):
				cell_id = data['cells'][i]
				center, _ = common.get_cell_geom(data, cell_id)
				data['centers'][cell_id] = center

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
