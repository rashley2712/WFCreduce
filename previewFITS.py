#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads a list of FITS files and previews them.')
	parser.add_argument('inputfiles', type=str, help='A text file listing which files to load.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-b', '--bias', type=str, help="Name of the bias frame.")
	parser.add_argument('-f', '--balance', type=str, help="Name of the balance (flat) frame.")
	parser.add_argument('-p', '--pause', type=float, default=0.05, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-j', '--json', type=str, help="Save to a JSON file (specify filename).")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True

	fileList = []

	listFile = open(arg.inputfiles, 'rt')
	for line in listFile:
		filename = line.strip()
		fileList.append(filename)
	listFile.close()

	print(fileList)

	if arg.bias is not None:
		print("loading the bias frame", arg.bias)
		hdul = astropy.io.fits.open(arg.bias)
		bias = hdul[1].data
		hdul.close()

	if arg.balance is not None:
		print("loading the balance frame", arg.balance)
		hdul = astropy.io.fits.open(arg.balance)
		balance = hdul[0].data
		hdul.close()

	
	for frameNo, FITSfile in enumerate(fileList):
		hdul = astropy.io.fits.open(FITSfile)
		#print(hdul.info())
		header = hdul[0].header
		for h in header:
			print(h, header[h])
		imageData = hdul[1].data
		# Subtract the bias
		if arg.bias is not None: imageData = imageData - bias
		# Divide by the balance frame (apply the flat)
		if arg.balance is not None: imageData = numpy.divide(imageData, balance)
		
		imageData = numpy.rot90(imageData)

		#print(numpy.shape(imageData))
		hdul.close()

		amplifiedImage = generallib.percentiles(imageData, 5, 95)
		matplotlib.pyplot.imshow(amplifiedImage)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.show(block=False)
		print("Frame number: {:d} run: {:s}".format(frameNo, FITSfile))
		matplotlib.pyplot.pause(arg.pause)
		matplotlib.pyplot.clf()
		#print(imageData)
	
	sys.exit()
