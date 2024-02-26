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

from os.path import basename, dirname, join
from sys import argv

DefineScalarExpression('number_density', 'mass_density / 1.67262192369e-21')
DefineScalarExpression('vx', 'velocity[0]')
DefineScalarExpression('vy', 'velocity[1]')
DefineScalarExpression('vz', 'velocity[2]')
DefineScalarExpression('Bx', 'perturbed_magnetic_field_cell_center[0]')
DefineScalarExpression('By', 'perturbed_magnetic_field_cell_center[1]')
DefineScalarExpression('Bz', 'perturbed_magnetic_field_cell_center[2]')

save_win_attrs = SaveWindowAttributes()
save_win_attrs.family = 0
save_win_attrs.format = save_win_attrs.PNG
save_win_attrs.resConstraint = save_win_attrs.NoConstraint

annot_attrs = AnnotationAttributes()
annot_attrs.userInfoFlag = 0
SetAnnotationAttributes(annot_attrs)


'''
Writes slice plot perpendicular to given axis of given variable.

Slice is taken at middle of grid.
Plot is saved in same dir as data.
'''
def save_slice_plot(outname, variable = 'mass_density', axis = 'x'):
	# get grid extents
	DeleteAllPlots()
	AddPlot('Mesh', 'mesh', 1, 1)
	DrawPlots()
	Query('SpatialExtents')
	xmin, xmax, ymin, ymax, zmin, zmax = GetQueryOutputValue()
	dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin

	# plot slice perpendicular to given axis of given variable
	DeleteAllPlots()
	AddPlot('Pseudocolor', variable, 1, 1)
	AddPlot('Mesh', 'mesh')
	AddOperator("Slice", 1)
	attrs = SliceAttributes()
	attrs.originType = attrs.Intercept
	attrs.project2d = 1
	if axis == 'x':
		attrs.originIntercept = dx / 2
		attrs.normal = (1, 0, 0)
		attrs.axisType = attrs.XAxis
		attrs.upAxis = (0, 0, 1)
		attrs.flip = 1
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dy/dz)
	elif axis == 'y':
		attrs.originIntercept = dy / 2
		attrs.normal = (0, 1, 0)
		attrs.axisType = attrs.YAxis
		attrs.upAxis = (0, 0, 1)
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dx/dz)
	elif axis == 'z':
		attrs.originIntercept = dz / 2
		attrs.normal = (0, 0, 1)
		attrs.axisType = attrs.ZAxis
		attrs.upAxis = (0, 1, 0)
		save_win_attrs.height = 900
		save_win_attrs.width = int(900*dx/dy)
	else:
		exit('Unsupported axis in save_slice_plot: ' + axis)
	SetOperatorOptions(attrs, 0, 1)
	save_win_attrs.fileName = outname
	SetSaveWindowAttributes(save_win_attrs)
	DrawPlots()
	SaveWindow()


for arg in argv[1:]:
	OpenDatabase(arg)
	for variable in [
		'mass_density',
		'number_density',
		'pressure',
		'vx', 'vy', 'vz',
		'Bx', 'By', 'Bz'
	]:
		for axis in ['x', 'y', 'z']:
			outname = join(dirname(arg), basename(arg).replace('.vtk', '_' + variable + '_' + axis + '.png'))
			save_slice_plot(outname, variable, axis)
exit()
