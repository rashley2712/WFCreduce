#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot, matplotlib.patches
import datetimelib
import photometrylib
import generallib
import astropy
import subprocess

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Makes a smaller flat for a windowed run.')
	parser.add_argument('flat', type=str, help='The name of the flat.')
	parser.add_argument('sample', type=str, help='A sample frame from the science run.')
	parser.add_argument("--preview",  action="store_true", help="Preview each CCD in 'ds9'.")
	parser.add_argument("-o", "--outputfilename", type=str, help="Output filename for the flat.")
	arg = parser.parse_args()

	
	hdul = astropy.io.fits.open(arg.flat)
	# Check if WFC
	numCCDs = 1
	header = hdul[0].header
	flatHeaders = {}
	flatCCDs = []
	for h in header:
		flatHeaders[h] = header[h]
	if len(hdul)==5:
		print("This exposure has multiple extensions... assuming WFC full image")
		numCCDs = len(hdul) - 1
	for CCD in range(1, numCCDs+1):
		flatCCDs.append(hdul[CCD].data)
	hdul.close()
	if numCCDs==4:
		flatData = numpy.array(flatCCDs[3])
	else: 
		flatData = numpy.array(flatCCDs[0])
	
	matplotlib.pyplot.figure()
	amplifiedImage = generallib.percentiles(flatData, 10, 98)
	matplotlib.pyplot.imshow(amplifiedImage)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.show(block=False)
	matplotlib.pyplot.pause(.05)

	# Load science frame
	hdul = astropy.io.fits.open(arg.sample)
	header = hdul[0].header
	sampleHeaders = {}
	for h in header:
		sampleHeaders[h] = header[h]
		#print(h, header[h])

	science = {"CCDXBIN": sampleHeaders['CCDXBIN'],
			   "CCDYBIN": sampleHeaders['CCDYBIN'],
			   "WINSEC4": sampleHeaders['WINSEC4']}

	print(science)
	science['data'] = hdul[1].data
	
	if 'enabled' in science['WINSEC4']:
		start = science['WINSEC4'].find('[')
		end = science['WINSEC4'].find(']')
		dimensionString = science['WINSEC4'][start+1:end]
		print(dimensionString)
		xDimensions = dimensionString.split(',')[0]
		xStart = int(xDimensions.split(':')[0])
		xEnd = int(xDimensions.split(':')[1])
		yDimensions = dimensionString.split(',')[1]
		yStart = int(yDimensions.split(':')[0])
		yEnd = int(yDimensions.split(':')[1])
		print(xStart, xEnd, yStart, yEnd)
		offset = 47
		offsety = -1
		sampledFlat = flatData[yStart+offsety:yEnd+offsety,xStart+offset:xEnd+offset+1]
		# Create a Rectangle patch
		rect = matplotlib.patches.Rectangle((xStart,yStart),(xEnd-xStart),(yEnd-yStart),linewidth=1,edgecolor='r',facecolor='none', ls=":")
		# Add the patch to the Axes
		
		matplotlib.pyplot.gca().add_patch(rect)
		matplotlib.pyplot.draw()
	print(numpy.shape(flatData))
	print(numpy.shape(sampledFlat))
	print(numpy.shape(science['data']))
	#sampledFlat = sampledFlat[::2, ::2]
	from skimage.measure import block_reduce
	sampledFlat = block_reduce(sampledFlat, block_size=(2, 2), func=numpy.mean)
	print(numpy.shape(sampledFlat))
	
	matplotlib.pyplot.figure()
	amplifiedImage = generallib.percentiles(sampledFlat, 10, 98)
	matplotlib.pyplot.imshow(amplifiedImage)
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.show(block=False)

	matplotlib.pyplot.figure()
	amplifiedImage = generallib.percentiles(science['data'], 10, 98)
	matplotlib.pyplot.imshow(amplifiedImage)
	matplotlib.pyplot.gca().invert_yaxis()
	
	deflatted = numpy.divide(science['data'],sampledFlat)
	matplotlib.pyplot.figure()
	amplifiedImage = generallib.percentiles(deflatted, 10, 98)
	matplotlib.pyplot.imshow(amplifiedImage)
	matplotlib.pyplot.gca().invert_yaxis()
	
	matplotlib.pyplot.show(block=False)
	matplotlib.pyplot.pause(.05)

	hdu = astropy.io.fits.PrimaryHDU(sampledFlat)
	hdul = astropy.io.fits.HDUList([hdu])
	hdul.writeto('sampled_balance.fits', overwrite=True)
		
	if arg.preview:
		ds9Command = ['ds9']
		ds9Command.append('sampled_balance.fits')
		ds9Command.append('-zscale')
		subprocess.Popen(ds9Command)

	input("Press enter to continue.")
	sys.exit()

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

			amplifiedImage = generallib.percentiles(imageData, 50, 95)
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
	else: flatFilename = "flat_%s.fits"%(filter)
	biasList.writeto(flatFilename, overwrite=True)
	print("Written flat to %s"%flatFilename)
	
	sys.exit()
