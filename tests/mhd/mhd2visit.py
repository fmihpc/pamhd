#! /usr/bin/env -S visit -nowin -cli -s
'''
Program for plotting MHD output of PAMHD using VisIt.

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

from argparse import ArgumentParser
from imp import load_source
from math import isnan
from os import _exit
from os.path import basename, dirname, join, realpath
from pathlib import Path
from sys import argv, stdout
common = load_source(
	'common',
	join(Path(realpath(__file__)).parent.parent, 'common.py')
)

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
DefineScalarExpression('B1x', 'volume_B1[0]')
DefineScalarExpression('B1y', 'volume_B1[1]')
DefineScalarExpression('B1z', 'volume_B1[2]')

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
	attrs.limitsMode = attrs.OriginalData
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
		save_win_attrs.width = int(min(2560, 900*dy/dz))
	elif dimension == 'y':
		attrs.originIntercept = ymin + dy / 2
		attrs.normal = (0, 1, 0)
		attrs.axisType = attrs.YAxis
		attrs.upAxis = (0, 0, 1)
		save_win_attrs.height = 900
		save_win_attrs.width = int(min(2560, 900*dx/dz))
	elif dimension == 'z':
		attrs.originIntercept = zmin + dz / 2
		attrs.normal = (0, 0, 1)
		attrs.axisType = attrs.ZAxis
		attrs.upAxis = (0, 1, 0)
		save_win_attrs.height = 900
		save_win_attrs.width = int(min(2560, 900*dx/dy))
	else:
		print('Unsupported dimension in save_slice_plot:', dimension)
		_exit(1)
	SetOperatorOptions(attrs, 0, 1)
	save_win_attrs.fileName = outname
	SetSaveWindowAttributes(save_win_attrs)
	DrawPlots()
	SaveWindow()
	attrs = PseudocolorAttributes()
	attrs.scaling = attrs.Linear
	SetPlotOptions(attrs)


plot_vars = {
	'mass_density', 'pressure', 'number_density', 'divB',
	'rank', 'V', 'Vx', 'Vy', 'Vz', 'B', 'Bx', 'By', 'Bz',
	'B0', 'B0x', 'B0y', 'B0z', 'B1', 'B1x', 'B1y', 'B1z',
	'target_ref_lvl_min', 'target_ref_lvl_max', 'Berror',
	'substep_period', 'cell_dt', 'mhd_info'
}
plot_vars_str = ''
for v in sorted(plot_vars):
	plot_vars_str += v + ','

parser = ArgumentParser()
parser.add_argument('files', metavar = 'F', nargs = '*', help = 'Names of input files in vtk format to plot')
parser.add_argument('--mesh', action = 'store_true', default = False, help = 'Plot mesh on top of variables (default: False)')
parser.add_argument('--dimensions', type = str, default = 'x,y,z', help = 'comma separated list of dimensions to plot (default: x,y,z)')
parser.add_argument('--variables', type = str, default = plot_vars_str, help = 'comma separated list of variables to plot (append _log for log plot, default: ' + plot_vars_str + ')')
parser.add_argument('--limits', type = str, default = '', help = 'comma separated list of range limits for variables given as var1_min,var1_max,,var2_max,,,var4_min...; with empty automatic limit(s), use / to escape - if var1_min < 0')
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

args.limits = args.limits.replace('/', '')
args.limits = args.limits.split(',')
while len(args.limits) < 2*len(args.variables):
	args.limits.append('')
for filename_ in args.files:
	if filename_.endswith('.dc'):
		with open(filename_, 'rb') as infile:
			data = common.get_metadata(infile)
			if data['file_version'] != 4:
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
		print("Couldn't load " + filename + ': ' + result)
		_exit(1)
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
_exit(0)
