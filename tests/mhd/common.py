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
	from numpy import dtype, fromfile, int32, uint64
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
	return ret_val


'''
Returns list of simulation cells and their byte offsets in infile.

Should be called after get_metadata() without changing infile read position.
'''
def get_cells(infile, total_cells):
	from numpy import dtype, fromfile
	return fromfile(infile, dtype = '2uint64', count = total_cells[0])


'''
Returns data of one cell from given byte offset in infile.
'''
def get_cell_data(infile, offset):
	from numpy import dtype, fromfile
	cell_data_t = 'double, 3double, double, 3double, intc, intc, 3double, 3double, 12double, 3double, 3double, 3double, 3double, 3double, 3double, double, intc'
	infile.seek(offset, 0)
	return fromfile(infile, dtype = cell_data_t, count = 1)
