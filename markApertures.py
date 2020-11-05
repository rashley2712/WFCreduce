#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy
import json
import classes
import shift


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Shows a list of apertures allows you to select the ones of interest.')
	parser.add_argument('inputfiles', type=str, help='A text file listing which files to process.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-b', '--bias', type=str, help="Name of the bias frame.")
	parser.add_argument('-f', '--balance', type=str, help="Name of the balance (flat) frame.")
	parser.add_argument('-p', '--pause', type=float, default=0.01, help="Number of seconds to pause on each plot.")
	parser.add_argument('-n', '--nframes', type=int, help="Number of frames to process before stopping. Default is all frames.")
	arg = parser.parse_args()
	if arg.nframes is not None:
		stopFrame = arg.nframes
	else:
		stopFrame = 1E8

	plotWidth = 8
	plotHeight = 8/1.7

	fileList = []

	listFile = open(arg.inputfiles, 'rt')
	for line in listFile:
		filename = line.strip()
		fileList.append(filename)
	listFile.close()

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


	rootApertures = classes.apertureDB()
	rootApertures.load()
	cat1 = numpy.array(rootApertures.makeCatalog())
	offsets = []
	frameList = classes.frameDB()
	for frameNo, FITSfile in enumerate(fileList):
		hdul = astropy.io.fits.open(FITSfile)
		FITSHeaders = {}
		header = hdul[0].header
		for h in header:
			FITSHeaders[h] = header[h]
		imageData = hdul[1].data
		hdul.close()
		# Subtract the bias
		if arg.bias is not None: imageData = imageData - bias
		# Divide by the balance frame (apply the flat)
		if arg.balance is not None: imageData = numpy.divide(imageData, balance)
		
		# Find the point sourcesBuild the average frame
		from astropy.stats import sigma_clipped_stats
		from photutils import datasets
		mean, median, std = sigma_clipped_stats(imageData, sigma=3.0)
		#print((mean, median, std))  
		from photutils import DAOStarFinder
		daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
		sources = daofind(imageData - median)  
		for col in sources.colnames:
			sources[col].info.format = '%.8g'  # for consistent table output
		#print(sources)
		apertureList = classes.apertureDB()
		cat2 = []
		for s in sources:
			target = classes.apertureClass((s['xcentroid'], s['ycentroid']), s['peak'])
			apertureList.add(target)
		apertureList.sort()
		apertureList.enableTrackers(n=10)
		cat2 = numpy.array(apertureList.makeCatalog())
		#print(cat2)
		
		dmax = 50
		fwhm = 4
		psize = 0.5
		mmax = 3
		img, xp,yp,xr,yr = shift.vimage(cat1, cat2, dmax, psize, fwhm)
		offset = { 'frame': frameNo, 'filesource': FITSfile, 'dx': xr, 'dy': yr}
		offsets.append(offset)
		
		# Draw the image
		amplifiedImage = numpy.rot90(generallib.percentiles(imageData, 5, 95))
		matplotlib.pyplot.imshow(amplifiedImage)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(arg.pause)
		matplotlib.pyplot.clf()
		print("Frame number: {:d}  ({:.1f}, {:.1f})   : {:s}".format(frameNo, offset['dx'], offset['dy'], FITSfile))
		frameInfo = { 'frame': frameNo, 'filesource': FITSfile, 'dx': xr, 'dy': yr}
		frameInfo['JD'] = FITSHeaders['JD']
		frameInfo['exptime'] = FITSHeaders['EXPTIME']
		frameInfo['JD'] = frameInfo['JD'] + frameInfo['exptime']/86400 / 2
		
		frameList.add(frameInfo)
		if frameNo>=stopFrame-1: break

	offsetPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	xValues = [o['frame'] for o in offsets]
	yValues = [o['dx'] for o in offsets]
	matplotlib.pyplot.scatter(xValues, yValues, label='x-offset')
	xValues = [o['frame'] for o in offsets]
	yValues = [o['dy'] for o in offsets]
	matplotlib.pyplot.scatter(xValues, yValues, label="y-offset")
	matplotlib.pyplot.legend()
	matplotlib.pyplot.xlabel("Frame number")
	matplotlib.pyplot.ylabel("Pixel shift")
	matplotlib.pyplot.plot([0, frameNo], [0, 0], ls=":")
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	
	frameList.save()

	input("Press enter to continue")
	sys.exit()
