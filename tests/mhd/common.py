#! /usr/bin/env python3

'''
Common functions for working with MHD output of PAMHD.

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


'''
Returns dict of simulation metadata in given file open for reading in binary mode.

Should be called before get_cells() without changing infile read position.
'''
def get_metadata(infile):
	from sys import version_info
	from numpy import fromfile
	ret_val = dict()
	file_version = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
	if file_version != 3:
		exit('Unsupported file version: ' + str(file_version))
	ret_val['file_version'] = file_version
	ret_val['sim_step'] = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
	ret_val['sim_time'], \
	ret_val['adiabatic_index'], \
	ret_val['proton_mass'], \
	ret_val['vacuum_permeability'] \
		= [float(i) for i in fromfile(infile, dtype = '4double', count = 1)[0]]
	nr_variables = int(fromfile(infile, dtype = 'u1', count = 1)[0])
	var_offsets = [int(i) for i in fromfile(infile, dtype = 'uint64', count = nr_variables)]
	endianness = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
	if endianness != 0x1234567890abcdef:
		exit('Unsupported endianness: ' + str(endianness))
	ret_val['endianness'] = endianness
	ret_val['ref_lvl_0_cells'] = [int(i) for i in fromfile(infile, dtype = '3uint64', count = 1)[0]]
	ret_val['max_ref_lvl'] = int(fromfile(infile, dtype = 'intc', count = 1)[0])
	ret_val['neighborhood_length'] = int(fromfile(infile, dtype = 'uintc', count = 1)[0])
	ret_val['periodicity'] = [int(i) for i in fromfile(infile, dtype = '3uint8', count = 1)[0]]
	geometry_id = int(fromfile(infile, dtype = 'intc', count = 1)[0])
	if geometry_id != 1:
		exit('Unsupported geometry: ' + str(geometry_id))
	ret_val['geometry_id'] = geometry_id
	ret_val['grid_start'] = [float(i) for i in fromfile(infile, dtype = '3double', count = 1)[0]]
	ret_val['lvl_0_cell_length'] = [float(i) for i in fromfile(infile, dtype = '3double', count = 1)[0]]
	ret_val['total_cells'] = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
	ret_val['cell_list_start'] = infile.tell()
	ret_val['var_data_start'] = dict()
	name_len = 8
	for offset in var_offsets:
		infile.seek(offset, 0)
		if version_info[0] >= 3:
			name = str(fromfile(infile, dtype = 'byte', count = name_len), 'ascii').strip()
		else:
			name = ''.join(fromfile(infile, dtype = 'c', count = name_len)).strip()
		ret_val['var_data_start'][name] = offset + name_len
	infile.seek(ret_val['cell_list_start'], 0)
	return ret_val


'''
Returns list of simulation cells and their byte offsets in infile.

Should be called after get_metadata() without changing infile read position.

Reads from current position of infile, to read fewer cells,
seek to place where previous read left off before calling again.
'''
def get_cells(infile, total_cells):
	from numpy import fromfile
	return [int(i) for i in fromfile(infile, dtype = 'uint64', count = total_cells)]


'''
Returns data of given range() of cells in cell list(s) returned by get_cells().

See end of file for example usage.
'''
def get_cell_data(infile, metadata, range_):
	from numpy import fromfile
	ret_val = []
	for i in range_:
		ret_val.append(dict())
		for varname in metadata['var_data_start']:
			infile.seek(metadata['var_data_start'][varname], 0)
			if varname == 'mhd':
				infile.seek(i*8*8, 1)
				ret_val[-1][varname] = fromfile(
					infile,
					dtype = 'double, 3double, double, 3double',
					count = 1)[0]
			elif varname == 'primary':
				infile.seek(i*(6+12), 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '18u1', count = 1)[0])
			elif varname == 'bgB':
				infile.seek(i*3*6*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '18double', count = 1)[0])
			elif varname == 'divfaceB':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'double', count = 1)[0]
			elif varname == 'edgeE':
				infile.seek(i*12*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '12double', count = 1)[0])
			elif varname == 'faceB':
				infile.seek(i*6*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '6double', count = 1)[0])
			elif varname == 'rank':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'mhd info':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'uintc', count = 1)[0]
			elif varname == 'ref lvls':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '2intc', count = 1)[0])
			elif varname == 'substep':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			else:
				if varname != 'fluxes':
					exit('Unsupported variable: ' + varname)
	return ret_val


def get_cell_geom(metadata, cell_id):
	orig_id = cell_id
	rl0c = metadata['ref_lvl_0_cells']
	l0cl = metadata['lvl_0_cell_length']
	gs = metadata['grid_start']
	mrl = metadata['max_ref_lvl']

	ref_lvl = 0
	while True:
		last = rl0c[0] * rl0c[1] * rl0c[2] * 2**(3*ref_lvl)
		if cell_id <= last:
			break
		cell_id -= last
		ref_lvl += 1
	if ref_lvl > mrl:
		raise Exception('Impossibly large refinement level for cell: ' + str(cell_id))

	cell_id -= 1
	cell_index = (
		cell_id % (rl0c[0] * 2**ref_lvl),
		(cell_id // (rl0c[0] * 2**ref_lvl)) % (rl0c[1] * 2**ref_lvl),
		cell_id // (rl0c[0]*rl0c[1] * 2**(2*ref_lvl)))
	cell_index = (
		cell_index[0] * 2**(mrl - ref_lvl),
		cell_index[1] * 2**(mrl - ref_lvl),
		cell_index[2] * 2**(mrl - ref_lvl))
	cell_center = ( # python2 division requires float(s)
		gs[0] + l0cl[0] * (cell_index[0] / 2.0**mrl + 1.0 / 2**(ref_lvl+1)),
		gs[1] + l0cl[1] * (cell_index[1] / 2.0**mrl + 1.0 / 2**(ref_lvl+1)),
		gs[2] + l0cl[2] * (cell_index[2] / 2.0**mrl + 1.0 / 2**(ref_lvl+1)))
	cell_length = (l0cl[0] / 2.0**ref_lvl, l0cl[1] / 2.0**ref_lvl, l0cl[2] / 2.0**ref_lvl)
	return cell_center, cell_length


if __name__ == '__main__':
	import sys
	for arg in sys.argv[1:]:
		infile = open(arg, 'rb')

		metadata = get_metadata(infile)
		print(metadata)

		total_cells = metadata['total_cells']
		cells1 = get_cells(infile, int(total_cells//2))
		pos = infile.tell()

		data1 = get_cell_data(infile, metadata, range(len(cells1)))
		if len(data1) > 0:
			print('Cell:', cells1[0], '\ndata:', data1[0])

		infile.seek(pos, 0)
		cells2 = get_cells(infile, int(total_cells - total_cells//2))

		data2 = get_cell_data(infile, metadata, range(len(cells1), len(cells1) + len(cells2)))
		if len(data2) > 0:
			print('Cell:', cells2[0], '\ndata:', data2[0])
