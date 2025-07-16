#! /usr/bin/env python3
'''
Common functions for working with MHD output of PAMHD.

Copyright 2024, 2025 Finnish Meteorological Institute

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
	if file_version != 4:
		exit('Unsupported file version: ' + str(file_version))
	ret_val['file_version'] = file_version
	ret_val['sim_step'] = int(fromfile(infile, dtype = 'uint64', count = 1)[0])
	ret_val['sim_time'], \
	ret_val['adiabatic_index'], \
	ret_val['proton_mass'], \
	ret_val['vacuum_permeability'], \
	ret_val['particle_temp_nrj_ratio'] \
		= [float(i) for i in fromfile(infile, dtype = '5double', count = 1)[0]]
	ret_val['nr_variables'] = int(fromfile(infile, dtype = 'u1', count = 1)[0])
	nr_variables = ret_val['nr_variables']
	ret_val['variable_offsets'] = [int(i) for i in fromfile(infile, dtype = 'uint64', count = nr_variables)]
	var_offsets = ret_val['variable_offsets']
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
			name = str(fromfile(infile, dtype = 'byte', count = name_len), 'ascii')
		else:
			name = ''.join(fromfile(infile, dtype = 'c', count = name_len))
		ret_val['var_data_start'][name] = offset + name_len
	infile.seek(ret_val['cell_list_start'], 0)
	return ret_val


'''
Writes simulation metadata to given file open for writing in binary mode.
'''
def set_metadata(outfile, metadata):
	from numpy import array
	if metadata['file_version'] != 3:
		raise Exception('Unsupported file version: ' + str(metadata['file_version']))
	if metadata['endianness'] != 0x1234567890abcdef:
		raise Exception('Unsupported endianness: ' + str(metadata['endianness']))
	if metadata['geometry_id'] != 1:
		raise Exception('Unsupported geometry: ' + str(metadata['geometry_id']))
	array(metadata['file_version'], dtype = 'uint64').tofile(outfile)
	array(metadata['sim_step'], dtype = 'uint64').tofile(outfile)
	array(metadata['sim_time'], dtype = 'double').tofile(outfile)
	array(metadata['adiabatic_index'], dtype = 'double').tofile(outfile)
	array(metadata['proton_mass'], dtype = 'double').tofile(outfile)
	array(metadata['vacuum_permeability'], dtype = 'double').tofile(outfile)
	array(metadata['nr_variables'], dtype = 'u1').tofile(outfile)
	array(metadata['variable_offsets'], dtype = 'uint64').tofile(outfile)
	array(metadata['endianness'], dtype = 'uint64').tofile(outfile)
	array(metadata['ref_lvl_0_cells'], dtype = 'uint64').tofile(outfile)
	array(metadata['max_ref_lvl'], dtype = 'intc').tofile(outfile)
	array(metadata['neighborhood_length'], dtype = 'uintc').tofile(outfile)
	array(metadata['periodicity'], dtype = 'uint8').tofile(outfile)
	array(metadata['geometry_id'], dtype = 'intc').tofile(outfile)
	array(metadata['grid_start'], dtype = 'double').tofile(outfile)
	array(metadata['lvl_0_cell_length'], dtype = 'double').tofile(outfile)
	array(metadata['total_cells'], dtype = 'uint64').tofile(outfile)


'''
Returns list of simulation cells in infile.

Should be called after get_metadata() without changing infile read position.

Reads from current position of infile, to read fewer cells at a time,
seek to place where previous read left off before calling again.
'''
def get_cells(infile, total_cells):
	from numpy import fromfile
	return [int(i) for i in fromfile(infile, dtype = 'uint64', count = total_cells)]


'''
Writes given list of simulation cells to outfile.

Should be called after set_metadata() without changing outfile write position.
'''
def set_cells(outfile, cells):
	from numpy import array
	array(cells, dtype = 'uint64').tofile(outfile)


'''
Returns data of given range() of cells in cell list(s) returned by get_cells().

See end of file for example usage.
'''
def get_cell_data(infile, metadata, range_, variables = None):
	from numpy import fromfile
	ret_val = []
	for i in range_:
		ret_val.append(dict())
		if variables == None:
			variables = list(metadata['var_data_start'])
		for varname in variables:
			infile.seek(metadata['var_data_start'][varname], 0)
			if varname == 'mhd     ':
				infile.seek(i*8*8, 1)
				ret_val[-1][varname] = fromfile(
					infile,
					dtype = 'double, 3double, double, 3double',
					count = 1)[0]
			elif varname == 'bgB     ':
				infile.seek(i*3*6*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '18double', count = 1)[0])
			elif varname == 'divfaceB':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'double', count = 1)[0]
			elif varname == 'faceB   ':
				infile.seek(i*6*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '6double', count = 1)[0])
			elif varname == 'rank    ':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'mhd info':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'ref lvls':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = list(fromfile(infile, dtype = '2intc', count = 1)[0])
			elif varname == 'substep ':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'substmin':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'substmax':
				infile.seek(i*4, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'intc', count = 1)[0]
			elif varname == 'timestep':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'double', count = 1)[0]
			elif varname == 'volE    ':
				infile.seek(i*3*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = '3double', count = 1)[0]
			elif varname == 'volJ    ':
				infile.seek(i*3*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = '3double', count = 1)[0]
			elif varname == 'nr ipart':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'uint64', count = 1)[0]
			elif varname == 'Berror  ':
				infile.seek(i*8, 1)
				ret_val[-1][varname] = fromfile(infile, dtype = 'double', count = 1)[0]
	return ret_val


'''
Writes simulation data to outfile.

Should be called after set_cells() without changing outfile write position.
Data should be in same order as cells given to set_cells().

Returns names and starting offsets of saved variables.
'''
def set_cell_data(outfile, metadata, data):
	from numpy import array

	ret_val = []
	if len(data) == 0:
		return ret_val

	for varname in metadata['var_data_start']:
		if varname == 'fluxes  ':
			continue

		ret_val.append((varname, outfile.tell()))
		outfile.write(varname.encode('utf-8'))
		if varname == 'mhd     ':
			for i in range(len(data)):
				try:
					d = data[i][varname]
				except Exception as e:
					print(varname, len(simdata), len(simcells), len(outcells), ind, e)
					exit(1)
				array((d[0], d[1][0], d[1][1], d[1][2], d[2], d[3][0], d[3][1], d[3][2]), dtype = 'double').tofile(outfile)
		if varname == 'bgB     ':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'double').tofile(outfile)
		if varname == 'divfaceB':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'double').tofile(outfile)
		if varname == 'edgeE   ':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'double').tofile(outfile)
		if varname == 'faceB   ':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'double').tofile(outfile)
		if varname == 'rank    ':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'intc').tofile(outfile)
		if varname == 'mhd info':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'uintc').tofile(outfile)
		if varname == 'ref lvls':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'intc').tofile(outfile)
		if varname == 'substep ':
			for i in range(len(data)):
				d = data[i][varname]
				array(d, dtype = 'intc').tofile(outfile)

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


# writes metadata in ascii vtk format to file open for writing
def write_metadata_vtk(outfile, meta):
	cells = meta['cells']
	outfile.write(
		"# vtk DataFile Version 2.0\n"
		+ "MHD data from PAMHD\n"
		+ "ASCII\nDATASET UNSTRUCTURED_GRID\n"
		+ "FIELD FieldData 2\n"
		+ "CYCLE 1 1 int\n" + str(meta['sim_step'])
		+ "\nTIME 1 1 double\n" +str(meta['sim_time'])
		+ "\nPOINTS " + str(8*len(cells)) + " float\n")
	if len(cells) == 0:
		return
	for c in cells:
		center, length = get_cell_geom(meta, c)
		start = (
			str(center[0] - length[0]/2) + ' ',
			str(center[1] - length[1]/2) + ' ',
			str(center[2] - length[2]/2) + ' ')
		end = (
			str(center[0] + length[0]/2) + ' ',
			str(center[1] + length[1]/2) + ' ',
			str(center[2] + length[2]/2) + ' ')
		outfile.write(
			start[0] + start[1] + start[2] + '\n'
			+ end[0] + start[1] + start[2] + '\n'
			+ start[0] + end[1] + start[2] + '\n'
			+ end[0] + end[1] + start[2] + '\n'
			+ start[0] + start[1] + end[2] + '\n'
			+ end[0] + start[1] + end[2] + '\n'
			+ start[0] + end[1] + end[2] + '\n'
			+ end[0] + end[1] + end[2] + '\n')

	outfile.write('CELLS ' + str(len(cells)) + ' ' + str(9*len(cells)) + '\n')
	for i in range(len(cells)):
		outfile.write('8 ')
		for j in range(8):
			outfile.write(str(i*8 + j) + ' ')
		outfile.write('\n')
	outfile.write('CELL_TYPES ' + str(len(cells)) + '\n')
	for i in range(len(cells)):
		outfile.write('11\n')
	outfile.write('CELL_DATA ' + str(len(cells)) + '\n')


# call after write_metadata_vtk
def write_variable_vtk(infile, outfile, variable, meta, nr_cells = -1):
	cells = meta['cells']
	chunks = []
	if nr_cells < 0 or nr_cells >= len(cells):
		chunks.append(range(len(cells)))
	else:
		start = 0
		end = nr_cells
		while start < len(cells):
			chunks.append(range(start, end))
			start = end
			end = min(end + nr_cells, len(cells))

	if variable == 'mhd     ':
		outfile.write('SCALARS mass_density double 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['mhd     '][0]) + '\n')
			del data
		outfile.write('VECTORS velocity double\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				mas = d['mhd     '][0]
				mom = d['mhd     '][1]
				outfile.write(
					str(mom[0]/mas) + ' '
					+ str(mom[1]/mas) + ' '
					+ str(mom[2]/mas) + '\n')
			del data
		outfile.write('SCALARS pressure double 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				mas = d['mhd     '][0]
				mom = d['mhd     '][1]
				nrj = d['mhd     '][2]
				mag = d['mhd     '][3]
				kin_nrj = (mom[0]**2 + mom[1]**2 + mom[2]**2) / 2 / mas
				mag_nrj = (mag[0]**2 + mag[1]**2 + mag[2]**2) / 2 / meta['vacuum_permeability']
				pressure = (nrj - kin_nrj - mag_nrj) * (meta['adiabatic_index'] - 1)
				outfile.write(str(pressure) + '\n')
			del data
		outfile.write('VECTORS volume_B1 double\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				mag = d['mhd     '][3]
				outfile.write(str(mag[0]) + ' ' + str(mag[1]) + ' ' + str(mag[2]) + '\n')
			del data

	if variable == 'bgB     ':
		outfile.write('VECTORS background_B double\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(
					str(d['bgB     '][3]) + ' '
					+ str(d['bgB     '][10]) + ' '
					+ str(d['bgB     '][17]) + '\n')
			del data

	if variable == 'totB    ':
		outfile.write('VECTORS total_magnetic_field double\n')
		for chunk in chunks:
			# write dont_solve cells as all 0
			data = get_cell_data(infile, meta, chunk, ['mhd     ', 'bgB     ', 'mhd info'])
			for d in data:
				if 'mhd info' in d and d['mhd info'] == 1:
					outfile.write('0 0 0\n')
					continue
				B0 = d['bgB     ']
				B1 = d['mhd     '][3]
				outfile.write(
					str(B1[0] + B0[3]) + ' '
					+ str(B1[1] + B0[10]) + ' '
					+ str(B1[2] + B0[17]) + '\n')
			del data

	if variable == 'divfaceB':
		outfile.write('SCALARS divergence_of_magnetic_field double 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['divfaceB']) + '\n')
			del data

	if variable == 'rank    ':
		outfile.write('SCALARS rank int 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['rank    ']) + '\n')
			del data

	if variable == 'ref lvls':
		outfile.write('SCALARS target_ref_lvl_min int 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['ref lvls'][0]) + '\n')
			del data
		outfile.write('SCALARS target_ref_lvl_max int 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['ref lvls'][1]) + '\n')
			del data

	if variable == 'substep ':
		outfile.write('SCALARS substep_period int 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['substep ']) + '\n')
			del data

	if variable == 'mhd info':
		outfile.write('SCALARS mhd_info int 1\nlookup_table default\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['mhd info']) + '\n')
			del data

	if variable == 'volE    ':
		outfile.write('VECTORS volume_E double\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(
					str(d['volE    '][0]) + ' '
					+ str(d['volE    '][1]) + ' '
					+ str(d['volE    '][2]) + '\n')
			del data

	if variable == 'volJ    ':
		outfile.write('VECTORS volume_J double\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(
					str(d['volJ    '][0]) + ' '
					+ str(d['volJ    '][1]) + ' '
					+ str(d['volJ    '][2]) + '\n')
			del data

	if variable == 'nr ipart':
		outfile.write('VECTORS nr_particles_internal int\n')
		for chunk in chunks:
			data = get_cell_data(infile, meta, chunk, [variable])
			for d in data:
				outfile.write(str(d['nr ipart']) + '\n')
			del data


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
