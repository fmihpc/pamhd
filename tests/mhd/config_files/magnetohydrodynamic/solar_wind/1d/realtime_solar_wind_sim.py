#! /usr/bin/env python
'''
Converts realtime ACE data from NOAA to format suitable for MHD test program of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names of the copyright holders nor the names of their contributors
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

from argparse import ArgumentParser
from datetime import datetime
from ftplib import FTP

def get_plasma_data():
	mag_lines = []
	swepam_lines = []
	conn = FTP('ftp.swpc.noaa.gov')
	conn.login()
	conn.retrlines('RETR pub/lists/ace/ace_mag_1m.txt', mag_lines.append)
	conn.retrlines('RETR pub/lists/ace/ace_swepam_1m.txt', swepam_lines.append)

	mag_data = []
	for line in mag_lines:
		line = line.strip()
		if line.startswith(':') or line.startswith('#'):
			continue

		yyyy, mm, dd, hhmm, mjd, sof, status, bx, by, bz, bt, lat, lon = line.split()

		if status != '0':
			continue

		if bx == '-999.9' or by == '-999.9' or bz == '-999.9':
			continue

		dt = datetime(int(yyyy), int(mm), int(dd), int(hhmm[0:2]), int(hhmm[2:4]))
		bx, by, bz = float(bx)*1e-9, float(by)*1e-9, float(bz)*1e-9
		mag_data.append((dt, bx, by, bz))

	# use average value of Bx
	Bx_avg = 0
	for item in mag_data:
		Bx_avg += item[1]
	Bx_avg /= len(mag_data)

	swepam_data = []
	for line in swepam_lines:
		line = line.strip()
		if line.startswith(':') or line.startswith('#'):
			continue

		yyyy, mm, dd, hhmm, mjd, sof, status, density, speed, temperature = line.split()

		if status != '0':
			continue

		if density == '-9999.9' or speed == '-9999.9' or temperature == '-1.00e+05':
			continue

		dt = datetime(int(yyyy), int(mm), int(dd), int(hhmm[0:2]), int(hhmm[2:4]))
		density, speed, temperature = float(density)*1e6, -float(speed)*1e3, float(temperature)
		if density <= 0:
			density = 5e4
		swepam_data.append((dt, density, speed, temperature ))

	merged_data = []
	mag_i = 0
	swepam_i = 0
	while mag_i < len(mag_data) and swepam_i < len(swepam_data):
		if mag_data[mag_i][0] == swepam_data[swepam_i][0]:
			merged_data.append((mag_data[mag_i][0:1] + (Bx_avg,) + mag_data[mag_i][2:] + swepam_data[swepam_i][1:]))
			mag_i += 1
			swepam_i += 1
		elif mag_data[mag_i][0] < swepam_data[swepam_i][0]:
			mag_i += 1
		elif mag_data[mag_i][0] > swepam_data[swepam_i][0]:
			swepam_i += 1

	return merged_data


def get_pressure(number_density, temperature):
	particle_temp_nrj_ratio = 1.3806488e-23
	return number_density * temperature * particle_temp_nrj_ratio


if __name__ == '__main__':

	parser = ArgumentParser(
		description
			= 'Convert real-time solar wind from NOAA into a 1d '
			+ 'run configuration for MHD test of PAMHD'
	)
	parser.add_argument(
		'--config',
		metavar = 'C',
		default = '',
		help = 'Write configuration to file with name C (if empty '
			+ 'use first time stamp of data in ISO 8601 format)'
	)
	parser.add_argument(
		'--model-output',
		metavar = 'O',
		default = '',
		help = 'Output model results into directory O (if empty '
			+ 'use first time stamp of data in near ISO 8601 format)'
	)
	parser.add_argument(
		'--cells',
		metavar = 'N',
		type = int,
		default = 1000,
		help = 'Use N simulation cells'
	)
	parser.add_argument(
		'--length',
		metavar = 'L',
		type = float,
		default = 1e9,
		help = 'Use simulation box of length L (m)'
	)
	parser.add_argument(
		'--duration',
		metavar = 'T',
		type = float,
		default = 1e4,
		help = 'Simulate solar wind for duration T (s)'
	)

	args = parser.parse_args()

	plasma_data = get_plasma_data()
	time_start = plasma_data[0][0]
	time_length = plasma_data[-1][0] - time_start
	time_length = time_length.days * 60*60*24 + time_length.seconds


	time_stamps_str = '\t\t\t"time-stamps": ['
	density_str = '\t\t\t"values": ['
	velocity_str = '\t\t\t"values": ['
	pressure_str = '\t\t\t"values": ['
	mag_str = '\t\t\t"values": ['
	for item in plasma_data:
		time = (item[0] - time_start).total_seconds()
		#time = time.days * 60*60*24 + time.seconds
		time_stamps_str += str(time) + ', '
		density_str += str(item[4]) + ', '
		velocity_str += '[' + str(item[5]) + ', 0, 0], '
		pressure_str += str(get_pressure(item[4], item[6])) + ', '
		mag_str += '[' + str(item[1]) + ', ' + str(item[2]) + ', ' + str(item[3]) + '], '
	time_stamps_str = time_stamps_str[:-2] + '],\n'
	density_str = density_str[:-2] + ']\n'
	velocity_str = velocity_str[:-2] + ']\n'
	pressure_str = pressure_str[:-2] + ']\n'
	mag_str = mag_str[:-2] + ']\n'

	config_file = None
	if args.config == '':
		config_file = open('config-' + time_start.isoformat() + '.json', 'w')
	else:
		config_file = open(args.config, 'w')

	if args.model_output == '':
		config_file.write('{\n"output-directory": "' + time_start.isoformat().replace(':', ''))
	else:
		config_file.write('{\n"output-directory": "' + args.model_output)

	config_file.write(
		'",\n"solver-mhd": "roe-athena",\n'
		+ '"time-start": 0,\n'
		+ '"time-length": ' + str(args.duration)
		+ ',\n"load-balancer": "RCB",\n'
		+ '"save-mhd-n": 60,\n'
		+ '"adiabatic-index": 1.6666666666666667,\n'
		+ '"vacuum-permeability": 1.2566370614359173e-06,\n'
		+ '"proton-mass": 1.6726217770000001e-27,\n'
		+ '"time-step-factor": 0.5,\n'
		+ '"grid-options": {\n'
		+ '\t"periodic": "{false, false, false}",\n'
		+ '\t"cells": "{100 + 2, 1, 1}",\n'
		+ '\t"volume": "{1e9 * (1 + 2 / (cells[0] - 2)), 1e9, 1e9}",\n'
		+ '\t"start": "{-1 * 1e9 / (cells[0] - 2), -volume[1]/2, -volume[2]/2}"\n'
		+ '},\n'
		+ '"geometries": [\n'
		+ '\t{"box": {\n'
		+ '\t\t"start": [-99e99, -99e99, -99e99],\n'
		+ '\t\t"end": [0, 99e99, 99e99]\n'
		+ '\t}},\n'
		+ '\t{"box": {\n'
		+ '\t\t"start": [1e9, -99e99, -99e99],\n'
		+ '\t\t"end": [99e99, 99e99, 99e99]\n'
		+ '\t}}\n],\n'
		+ '"number-density": {\n'
		+ '\t"default": ' + str(plasma_data[0][4]) + ',\n'
		+ '\t"copy-boundaries": [{"geometry-id": 0}],\n'
		+ '\t"value-boundaries": [\n\t\t{\n\t\t\t"geometry-id": 1,\n'
		+ time_stamps_str + density_str
		+ '\t\t}\n\t]\n},\n'
		+ '"velocity": {\n'
		+ '\t"default": [' + str(plasma_data[0][5]) + ', 0, 0],\n'
		+ '\t"copy-boundaries": [{"geometry-id": 0}],\n'
		+ '\t"value-boundaries": [\n\t\t{\n\t\t\t"geometry-id": 1,\n'
		+ time_stamps_str + velocity_str
		+ '\t\t}\n\t]\n},\n'
		+ '"pressure": {\n'
		+ '\t"default": ' + str(get_pressure(item[4], item[6])) + ',\n'
		+ '\t"copy-boundaries": [{"geometry-id": 0}],\n'
		+ '\t"value-boundaries": [\n\t\t{\n\t\t\t"geometry-id": 1,\n'
		+ time_stamps_str + pressure_str
		+ '\t\t}\n\t]\n},\n'
		+ '"magnetic-field": {\n'
		+ '\t"default": [' + str(plasma_data[0][1]) + ', ' + str(plasma_data[0][2]) + ', ' + str(plasma_data[0][3]) + '],\n'
		+ '\t"copy-boundaries": [{"geometry-id": 0}],\n'
		+ '\t"value-boundaries": [\n\t\t{\n\t\t\t"geometry-id": 1,\n'
		+ time_stamps_str + mag_str
		+ '\t\t}\n\t]\n}}\n'
	)
