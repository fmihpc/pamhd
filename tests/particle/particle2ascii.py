#! /usr/bin/env python3
'''
Converts particle data of PAMHD to ASCII format.

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


def convert(inname, verbose):
	try:
		from numpy import fromfile
	except Exception as e:
		print("Couldn't import numpy.fromfile, try: pip install --user numpy: ", e)
		exit(1)

	try:
		infile = open(inname, 'rb')
	except Exception as e:
		raise Exception("Couldn't open file " + inname + ' for reading: ' + str(e))

	if verbose:
		print('Converting file', inname)

	metadata = common.get_metadata(infile)
	if metadata['file_version'] != 4:
		exit('Unsupported file version: ' + str(metadata['file_version']))
	if verbose:
		print('Simulation step:', metadata['sim_step'])
		print('Simulation time:', metadata['sim_time'])

	if metadata['geometry_id'] != 1:
		raise Exception('Unsupported geometry in file ' + inname)

	if verbose:
		print('Number of simulation cells:', metadata['total_cells'])
		print('Number of variables:', len(metadata['var_data_start']), '\nVariables:')
		for name in metadata['var_data_start']:
			print(name, 'start in file:', metadata['var_data_start'])

	metadata['cells'] = common.get_cells(infile, metadata['total_cells'])

	outname = inname.replace('.dc', '_part.txt')
	outfile = open(outname, 'w')

	if not 'nr ipart' in metadata['var_data_start']:
		raise Exception('Data for number of particles not found in file ' + inname)
	if not 'ipart   ' in metadata['var_data_start']:
		raise Exception('Particle data not found in file ' + inname)
	# track byte offset of cell's particle data
	offset = metadata['var_data_start']['ipart   ']
	for i in range(len(metadata['cells'])):
		infile.seek(i*8 + metadata['var_data_start']['nr ipart'], 0)
		nr_ipart = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
		infile.seek(offset, 0)
		# position, velocity, mass, species mass, charge/mass, id
		data = fromfile(
			infile,
			dtype = '3double, 3double, double, double, double, uint64',
			count = nr_ipart)
		for d in data:
			pos, vel, mas, spm, c2m, id_ = d
			outfile.write(str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2]) + ' ')
			outfile.write(str(vel[0]) + ' ' + str(vel[1]) + ' ' + str(vel[2]) + ' ')
			outfile.write(str(mas) + ' ' + str(spm) + ' ' + str(c2m) + ' ' + str(id_) + '\n')
		offset += nr_ipart*8*(3+3+1+1+1+1)

if __name__ == '__main__':
	import sys
	verbose = False
	for arg in sys.argv[1:]:
		if arg == '--verbose':
			verbose = True
			break
	for arg in sys.argv[1:]:
		convert(arg, verbose)
