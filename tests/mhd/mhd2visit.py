#! /usr/bin/env -S visit -nowin -cli -s
'''
Program for plotting MHD output of PAMHD using VisIt.

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

DefineScalarExpression('number_density', 'mass_density / 1.67262192369e-21')
DefineScalarExpression('V', 'velocity_magnitude')
DefineScalarExpression('Vx', 'velocity[0]')
DefineScalarExpression('Vy', 'velocity[1]')
DefineScalarExpression('Vz', 'velocity[2]')
DefineScalarExpression('divB', 'divergence_of_magnetic_field')
DefineScalarExpression('B', 'total_face_magnetic_field_magnitude')
DefineScalarExpression('Bx', 'total_face_magnetic_field[0]')
DefineScalarExpression('By', 'total_face_magnetic_field[1]')
DefineScalarExpression('Bz', 'total_face_magnetic_field[2]')
DefineScalarExpression('B0', 'background_B_magnitude')
DefineScalarExpression('B0x', 'background_B[0]')
DefineScalarExpression('B0y', 'background_B[1]')
DefineScalarExpression('B0z', 'background_B[2]')
DefineScalarExpression('B1', 'volume_B_magnitude')
DefineScalarExpression('B1x', 'volume_B[0]')
DefineScalarExpression('B1y', 'volume_B[1]')
DefineScalarExpression('B1z', 'volume_B[2]')

save_win_attrs = SaveWindowAttributes()
save_win_attrs.family = 0
save_win_attrs.format = save_win_attrs.PNG
save_win_attrs.resConstraint = save_win_attrs.NoConstraint

annot_attrs = AnnotationAttributes()
annot_attrs.userInfoFlag = 0
SetAnnotationAttributes(annot_attrs)


'''
Writes slice plot of variable perpendicular to dimension.

Slice is taken at middle of grid.
Plot is saved in same dir as data.
'''
def save_slice_plot(
	outname,
	variable = 'mass_density_log',
	dimension = 'x',
	mesh = False,
	var_min = None,
	var_max = None
):
	# get grid extents
	DeleteAllPlots()
	AddPlot('Mesh', 'mesh', 1, 1)
	DrawPlots()
	Query('SpatialExtents')
	xmin, xmax, ymin, ymax, zmin, zmax = GetQueryOutputValue()
	dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin

	# plot slice perpendicular to given dimension of given variable
	DeleteAllPlots()
	AddPlot('Pseudocolor', variable.replace('_log', ''), 1, 1)
	attrs = PseudocolorAttributes()
	attrs.limitsMode = attrs.CurrentPlot
	if variable.endswith('log'):
		attrs.scaling = attrs.Log
	if var_min != None:
		attrs.minFlag = 1
		attrs.min = var_min
	if var_max != None:
		attrs.maxFlag = 1
		attrs.max = var_max
	SetPlotOptions(attrs)
	if mesh:
		AddPlot('Mesh', 'mesh')
	AddOperator("Slice", 1)
	attrs = SliceAttributes()
	attrs.originType = attrs.Intercept
	attrs.project2d = 1
	if dimension == 'x':
		attrs.originIntercept = xmin + dx / 2
		attrs.normal = (1, 0, 0)
		attrs.axisType = attrs.XAxis
		attrs.upAxis = (0, 0, 1)
		attrs.flip = 1
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dy/dz)
	elif dimension == 'y':
		attrs.originIntercept = ymin + dy / 2
		attrs.normal = (0, 1, 0)
		attrs.axisType = attrs.YAxis
		attrs.upAxis = (0, 0, 1)
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dx/dz)
	elif dimension == 'z':
		attrs.originIntercept = zmin + dz / 2
		attrs.normal = (0, 0, 1)
		attrs.axisType = attrs.ZAxis
		attrs.upAxis = (0, 1, 0)
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dx/dy)
	else:
		exit('Unsupported dimension in save_slice_plot: ' + dimension)
	SetOperatorOptions(attrs, 0, 1)
	save_win_attrs.fileName = outname
	SetSaveWindowAttributes(save_win_attrs)
	DrawPlots()
	SaveWindow()
	attrs = PseudocolorAttributes()
	attrs.scaling = attrs.Linear
	SetPlotOptions(attrs)


def dc2vtk(outname, data):
	with open(outname, 'w') as outfile:
		cells = data['cells']
		outfile.write(
			"# vtk DataFile Version 2.0\n"
			+ "MHD data from PAMHD\n"
			+ "ASCII\nDATASET UNSTRUCTURED_GRID\n"
			+ "FIELD FieldData 2\n"
			+ "CYCLE 1 1 int\n" + str(data['sim_step'])
			+ "\nTIME 1 1 double\n" +str(data['sim_time'])
			+ "\nPOINTS " + str(8*len(cells)) + " float\n")
		if len(cells) == 0:
			return
		for c in cells:
			center, length = common.get_cell_geom(data, c)
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

		if 'mhd' in data[cells[0]]:
			outfile.write('SCALARS mass_density double 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['mhd'][0]) + '\n')
			outfile.write('VECTORS velocity double\n')
			for c in cells:
				mas = data[c]['mhd'][0]
				mom = data[c]['mhd'][1]
				outfile.write(
					str(mom[0]/mas) + ' '
					+ str(mom[1]/mas) + ' '
					+ str(mom[2]/mas) + '\n')
			outfile.write('SCALARS pressure double 1\nlookup_table default\n')
			for c in cells:
				mas = data[c]['mhd'][0]
				mom = data[c]['mhd'][1]
				nrj = data[c]['mhd'][2]
				mag = data[c]['mhd'][3]
				kin_nrj = (mom[0]**2 + mom[1]**2 + mom[2]**2) / 2 / mas
				mag_nrj = (mag[0]**2 + mag[1]**2 + mag[2]**2) / 2 / data['vacuum_permeability']
				pressure = (nrj - kin_nrj - mag_nrj) * (data['adiabatic_index'] - 1)
				outfile.write(str(pressure) + '\n')
			outfile.write('VECTORS volume_B double\n')
			for c in cells:
				mag = data[c]['mhd'][3]
				outfile.write(str(mag[0]) + ' ' + str(mag[1]) + ' ' + str(mag[2]) + '\n')

		if 'primary' in data[cells[0]]:
			outfile.write('SCALARS nr_primary_faces int 1\nlookup_table default\n')
			for c in cells:
				nr = 0
				for i in range(6):
					if data[c]['primary'][i] > 0:
						nr += 1
				outfile.write(str(nr) + '\n')
			outfile.write('SCALARS nr_primary_edges int 1\nlookup_table default\n')
			for c in cells:
				nr = 0
				for i in range(6, 18):
					if data[c]['primary'][i] > 0:
						nr += 1
				outfile.write(str(nr) + '\n')

		if 'bgB' in data[cells[0]]:
			outfile.write('VECTORS background_B double\n')
			for c in cells:
				outfile.write(
					str(data[c]['bgB'][3]) + ' '
					+ str(data[c]['bgB'][10]) + ' '
					+ str(data[c]['bgB'][17]) + '\n')

		if 'mhd' in data[cells[0]] and 'bgB' in data[cells[0]]:
			outfile.write('VECTORS total_face_magnetic_field double\n')
			for c in cells:
				B0 = data[c]['bgB']
				B1 = data[c]['mhd'][3]
				outfile.write(
					str(B1[0] + B0[3]) + ' '
					+ str(B1[1] + B0[10]) + ' '
					+ str(B1[2] + B0[17]) + '\n')

		if 'divfaceB' in data[cells[0]]:
			outfile.write('SCALARS divergence_of_magnetic_field double 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['divfaceB']) + '\n')

		if 'rank' in data[cells[0]]:
			outfile.write('SCALARS rank int 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['rank']) + '\n')

		if 'ref lvls' in data[cells[0]]:
			outfile.write('SCALARS target_ref_lvl_min int 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['ref lvls'][0]) + '\n')

		if 'ref lvls' in data[cells[0]]:
			outfile.write('SCALARS target_ref_lvl_max int 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['ref lvls'][1]) + '\n')

		if 'substep' in data[cells[0]]:
			outfile.write('SCALARS substep_period int 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['substep']) + '\n')

		if 'mhd info' in data[cells[0]]:
			outfile.write('SCALARS mhd_info int 1\nlookup_table default\n')
			for c in cells:
				outfile.write(str(data[c]['mhd info']) + '\n')


plot_vars = {
	'mass_density', 'pressure',
	'number_density', 'divB', 'rank',
	'V', 'Vx', 'Vy', 'Vz',
	'B', 'Bx', 'By', 'Bz',
	'B0', 'B0x', 'B0y', 'B0z',
	'B1', 'B1x', 'B1y', 'B1z',
	'target_ref_lvl_min', 'target_ref_lvl_max',
	'substep_period', 'mhd_info'
}
plot_vars_str = ''
for v in sorted(plot_vars):
	plot_vars_str += v + ','

parser = ArgumentParser()
parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of input files in vtk format to plot')
parser.add_argument('--mesh', action = 'store_true', default = False, help = 'Plot mesh on top of variables (default: False)')
parser.add_argument('--dimensions', type = str, default = 'x,y,z', help = 'comma separated list of dimensions to plot (default: x,y,z)')
parser.add_argument('--variables', type = str, default = plot_vars_str, help = 'comma separated list of variables to plot (append _log for log plot, default: ' + plot_vars_str + ')')
parser.add_argument('--limits', type = str, default = '', help = 'comma separated list of range limits for variables given as var1_min,var1_max,,var2_max,,,var4_min...; with empty automatic limit(s)')
parser.add_argument('--verbose', action = 'store_true', default = False)
args = parser.parse_args()

if args.verbose:
	stdout.write('Plotting files: ' + str(args.files) + '\n')

args.dimensions = [d for d in args.dimensions.split(',') if d in {'x','y','z'}]
if args.verbose:
	stdout.write('Plotting dimensions: ' + str(args.dimensions) + '\n')

args.variables = [v for v in args.variables.split(',') if v.replace('_log','') in plot_vars]
if args.verbose:
	stdout.write('Plotting variables: ' + str(args.variables) + '\n')

args.limits = args.limits.split(',')
while len(args.limits) < 2*len(args.variables):
	args.limits.append('')
for filename_ in args.files:
	if filename_.endswith('.dc'):
		with open(filename_, 'rb') as infile:
			data = common.get_metadata(infile)
			if data['file_version'] != 3:
				print('Unsupported file version:', data['file_version'])
				continue
			if data['geometry_id'] != 1:
				print('Unsupported geometry:', data['geometry_id'])
				continue
			data['cells'] = common.get_cells(infile, data['total_cells'])
			sim_data_ = common.get_cell_data(infile, data, range(len(data['cells'])))
			for i in range(len(data['cells'])):
				data[data['cells'][i]] = sim_data_[i]
		outname = filename_[:-2] + 'vtk'
		dc2vtk(outname, data)
		filename = outname
	else:
		filename = filename_
	OpenDatabase(filename)
	result = GetLastError()
	if result != '':
		exit("Couldn't load " + filename + ': ' + result)
	for var_i in range(len(args.variables)):
		var = args.variables[var_i]
		var_min = None
		if args.limits[2*var_i] != '':
			var_min = float(args.limits[2*var_i])
		var_max = None
		if args.limits[2*var_i+1] != '':
			var_max = float(args.limits[2*var_i+1])
		for dim in args.dimensions:
			outname = join(dirname(filename), basename(filename).replace('.vtk', '_' + var + '_' + dim + '.png'))
			save_slice_plot(outname, var, dim, mesh = args.mesh, var_min = var_min, var_max = var_max)
	DeleteAllPlots()
	CloseDatabase(filename)
exit()
