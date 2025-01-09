#! /usr/bin/env python3
'''
Prints basic information about test program output files.

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
from sys import argv, modules
spec = spec_from_file_location(
	'common',
	join(Path(realpath(__file__)).parent.parent, 'common.py')
)
common = module_from_spec(spec)
modules['common'] = common
spec.loader.exec_module(common)

if __name__ == '__main__':
	for arg in argv[1:]:
		infile = open(arg, 'rb')
		meta = common.get_metadata(infile)
		print('Output file:', arg)
		print('File version:', meta['file_version'])
		print('Simulation step:', meta['sim_step'])
		print('Simulation time:', meta['sim_time'])
		print('Adiabatic index:', meta['adiabatic_index'])
		print('Proton mass:', meta['proton_mass'])
		print('Vacuum permeability:', meta['vacuum_permeability'])
		print('Number of saved simulation variables:', meta['nr_variables'])
		print('Variable byte offsets: ', end = '')
		for var in meta['variable_offsets']:
			print(var, end = ', ')
		print('\nEndianness:', hex(meta['endianness']))
		rl0c = meta['ref_lvl_0_cells']
		print('Number of ref lvl 0 cells in each dim:', rl0c[0], rl0c[1], rl0c[2])
		print('Maximum refinement level of grid:', meta['max_ref_lvl'])
		print("Length of cells' neighborhood:", meta['neighborhood_length'])
		per = meta['periodicity']
		print('Periodicity:', per[0], per[1], per[2])
		print('Geometry id:', meta['geometry_id'])
		start = meta['grid_start']
		print('Starting coordinates of grid:', start[0], start[1], start[2])
		dr = meta['lvl_0_cell_length']
		print('Length of ref lvl 0 cells:', dr[0], dr[1], dr[2])
		print('Total number of saved cells:', meta['total_cells'])
		print('Names of saved variables: ', end = '')
		for var in meta['var_data_start']:
			print(var, end = ', ')
		print()
