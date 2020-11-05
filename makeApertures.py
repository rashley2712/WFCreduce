#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy
import json
import classes

def onclick(event):
	global ix, iy
	ix, iy = event.xdata, event.ydata
	return None

def labelAperture(number):
	matplotlib.pyplot.figure(fig.number)
	aperture = apertureList.getTargetByNumber(number)
	if aperture is not None: 
		print("Labeling: ", aperture)
		matplotlib.pyplot.text(aperture.rootPosition[0]+10, aperture.rootPosition[1], str(number), color="red", size="x-large")
	matplotlib.pyplot.draw()

def keypress(event):
	print('press', event.key)
	ix, iy = event.xdata, event.ydata
	print('x = %d, y = %d'%(ix, iy))
	try:
		keynumber = int(event.key)
		apertureList.match(ix, iy, keynumber)
		labelAperture(keynumber)
	except Exception as e:
		print("Press a number to assign it to an aperture.")
		print(e)
	sys.stdout.flush()
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads a list of FITS files and previews them.')
	parser.add_argument('inputfiles', type=str, help='A text file listing which files to load.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-n', '--nframes', type=int, default=5, help="Number of frames to average before finding targets. Default value: 5.")
	parser.add_argument('-b', '--bias', type=str, help="Name of the bias frame.")
	parser.add_argument('-f', '--balance', type=str, help="Name of the balance (flat) frame.")
	parser.add_argument('-p', '--pause', type=float, default=0.01, help="Number of seconds to pause on each plot.")
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
		balance = hdul[1].data
		hdul.close()

	medianFrameStack = []	
	for frameNo, FITSfile in enumerate(fileList):
		hdul = astropy.io.fits.open(FITSfile)
		# print(hdul.info())
		# Put all of the headers into a dictionary object
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
		
		# Build the average frame
		if frameNo==0: average = imageData
		else: average = numpy.add(average, imageData)
		medianFrameStack.append(imageData)
		
	
		amplifiedImage = generallib.percentiles(imageData, 5, 95)
		matplotlib.pyplot.imshow(amplifiedImage)
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(arg.pause)
		matplotlib.pyplot.clf()
		print("Frame number: {:d}".format(frameNo))
		if frameNo==arg.nframes-1: break

	
	average = numpy.divide(average, frameNo)
	medianFrame = numpy.median(medianFrameStack, axis=0)
	amplifiedImage = generallib.percentiles(medianFrame, 5, 95)
	
	fig = matplotlib.pyplot.figure()
	matplotlib.pyplot.imshow(amplifiedImage)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.draw()

	from astropy.stats import sigma_clipped_stats
	from photutils import datasets
	mean, median, std = sigma_clipped_stats(medianFrame, sigma=3.0)
	print((mean, median, std))  
	from photutils import DAOStarFinder
	daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
	sources = daofind(medianFrame - median)  
	for col in sources.colnames:
		sources[col].info.format = '%.8g'  # for consistent table output
	print(sources)
	apertureList = classes.apertureDB()
	for s in sources:
		target = classes.apertureClass((s['xcentroid'], s['ycentroid']), s['peak'])
		apertureList.add(target)
	from photutils import CircularAperture
	positions = numpy.transpose((sources['xcentroid'], sources['ycentroid']))

	apertures = CircularAperture(positions, r=10.)
	apertures.plot(color='red', lw=1.5, alpha=0.7)
	
	apertureList.sort()
	apertureList.enableTrackers(10)
	apertureList.save()


	mouseclick = fig.canvas.mpl_connect('button_press_event', onclick)
	keypress = fig.canvas.mpl_connect('key_press_event', keypress)	
	matplotlib.pyplot.show(block=True)

	print("Saving apertures")
	apertureList.save()

	sys.exit()
