#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy
import subprocess

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads a list of FITS files and previews them.')
	parser.add_argument('inputfiles', type=str, help='A text file listing which files to load.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-b', '--bias', type=str, help="Name of the bias frame.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument("--preview",  action="store_true", help="Preview each CCD in 'ds9'.")
	parser.add_argument("-o", "--outputfilename", type=str, help="Output filename for the flat.")
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
	

	# Check if WFC
	numCCDs = 1
	checkFITSFile = fileList[0]
	hdul = astropy.io.fits.open(checkFITSFile)
	print(hdul.info())
	header = hdul[0].header
	allHeaders = {}
	for h in header:
		print(h, header[h])
		allHeaders[h] = header[h]
	if len(hdul)==5:
		print("This exposure has multiple extensions... assuming WFC full image")
		numCCDs = len(hdul) - 1


	if arg.bias is not None:
		print("loading the bias frame", arg.bias)
		hdul = astropy.io.fits.open(arg.bias)
		numBiasCCDs = len(hdul) - 1
		if numBiasCCDs!=numCCDs:
			print("There is a mismatch in the number of FITS extensions between the bias and the flat frames.")
			sys.exit()
		biases = []
		for i in range(1, numBiasCCDs+1):
			biases.append(hdul[i].data)
		hdul.close()


	for CCD in range(1, numCCDs+1):
		flatFrames = []
		for index, FITSfile in enumerate(fileList):
			frameNo = index+1
			print("CCD %d  Frame no: %d"%(CCD, frameNo))
			hdul = astropy.io.fits.open(FITSfile)
			header = hdul[0].header
			imageData = hdul[CCD].data
			# Subtract the bias
			if arg.bias is not None: imageData = imageData - biases[CCD-1]
			flatFrames.append(imageData)
			
			hdul.close()

			amplifiedImage = generallib.percentiles(imageData, 5, 95)
			matplotlib.pyplot.imshow(amplifiedImage)
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.show(block=False)
			matplotlib.pyplot.pause(.05)
			matplotlib.pyplot.clf()
		
		flat = numpy.median(flatFrames, axis = 0)
		amplifiedImage = generallib.percentiles(flat, 5, 95)
		matplotlib.pyplot.imshow(amplifiedImage)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.show(block=False)

		hdu = astropy.io.fits.PrimaryHDU(flat)
		hdul = astropy.io.fits.HDUList([hdu])
		hdul.writeto('flat_%d.fits'%CCD, overwrite=True)
		midx, midy = ( int(numpy.shape(flat)[0] / 2), int(numpy.shape(flat)[1] / 2))
		pixelRange = 15
		centralRegion = flat[midx-pixelRange:midx+pixelRange, midy-pixelRange:midy+pixelRange]
		print("Size of central sample:", numpy.shape(centralRegion),"or", numpy.shape(centralRegion)[0] * numpy.shape(centralRegion)[1],"pixels.")
		mean = numpy.mean(centralRegion)
		balance = numpy.divide(flat, mean)
		hdu = astropy.io.fits.PrimaryHDU(balance)
		hdul = astropy.io.fits.HDUList([hdu])
		hdul.writeto('balance_%d.fits'%CCD, overwrite=True)
		
		if arg.preview:
			ds9Command = ['ds9']
			ds9Command.append('flat_%d.fits'%CCD)
			ds9Command.append('-zscale')
			subprocess.Popen(ds9Command)

		#ds9Command = ['ds9']
		#ds9Command.append('balance.fits')
		#ds9Command.append('-zscale')
		#subprocess.Popen(ds9Command)
			

	# Reformat all CCDs into 1 fits file
	print("Compiling CCDs into a single balance .fits file.")
	hdr = astropy.io.fits.Header()
	masterhdu = astropy.io.fits.PrimaryHDU(header=hdr)
	biasList =  astropy.io.fits.HDUList([masterhdu])
	for CCD in range(1, numCCDs+1):
		print("CCD", CCD)
		hdul = astropy.io.fits.open('balance_%d.fits'%CCD)
		imageData = hdul[0].data
		biasList.append(astropy.io.fits.ImageHDU(imageData))
	filter = 'g'
	if arg.outputfilename is not None: flatFilename = arg.outputfilename
	else: flatFilename = "balance_%s.fits"%(filter)
	biasList.writeto(flatFilename, overwrite=True)
	print("Written flat to %s"%flatFilename)
	
	sys.exit()
