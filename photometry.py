#!/usr/bin/env python3
import argparse, sys, numpy
import scipy.optimize
import matplotlib.pyplot
import datetimelib
import photometrylib
import generallib
import astropy
import json
import classes
import shift


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Compute photometry across all frames')
	parser.add_argument('--framedb', type=str, default="frames.json", help='The JSON file containing the frame list.')
	parser.add_argument('--apertures', type=str, default="apertures.json", help="The JSON file containing the apertures.")
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



	frameList = classes.frameDB()
	frameList.load()
	apertures = classes.apertureDB()
	apertures.load()
	targets = apertures.getTargets()
	targetList = []
	for index, t in enumerate(targets):
		print(index, t)
		target = classes.target(id=index)
		targetList.append(target)
	print("Number of reduction apertures:", len(targets))
	
	plot = False
	if plot: 
		cutFigure = matplotlib.pyplot.figure()
		fitFigure = matplotlib.pyplot.figure()
		plotGaussian = matplotlib.pyplot.figure()
	for frame in frameList.allFrames:
		FITSfile = frame['filesource']
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
		

		positions = ([t.rootPosition[0] + frame['dx'] for t in targets], [t.rootPosition[1] + frame['dy'] for t in targets])
		positions = numpy.transpose(positions)
		
		from photutils import CircularAperture
		from photutils import CircularAnnulus
		from photutils import aperture_photometry
		from photutils.utils import calc_total_error

		# Fit a Gaussian 2-D to the first aperture
		# Cut out a bit of the image near the target
		t = targets[0]
		size = 20
		offsety = -int(frame['dy'])
		offsetx = -int(frame['dx'])

		left, right = (int(t.rootPosition[0]-size/2), int(t.rootPosition[0]+size/2))
		top, bottom= (int(t.rootPosition[1]-size/2), int(t.rootPosition[1]+size/2))
		imageCut = imageData[ top-offsety:bottom-offsety , left-offsetx:right-offsetx ] 
		if plot: 
			from mpl_toolkits.mplot3d import Axes3D 
			from matplotlib import cm
			matplotlib.pyplot.figure(cutFigure.number)
			amplifiedImage = generallib.percentiles(imageCut, 5, 95)
			matplotlib.pyplot.imshow(amplifiedImage)
			matplotlib.pyplot.gca().invert_yaxis()
		x = numpy.arange(0, size, 1)
		y = numpy.arange(0, size, 1)
		x, y = numpy.meshgrid(x, y)
		z = imageCut
		
		def Gaussian(xy, amp, x0, y0, c):
			x, y = xy
			inner = (x - x0)**2 
			inner += (y - y0)**2
			#print(amp, x0, y0, c)
			return amp * numpy.exp(-inner/(2*c**2))

		FWHM = 10
		guess = [900, size/2, size/2, 4]
		xydata = numpy.vstack((x.ravel(),y.ravel()))
		zdata = z.ravel()
		popt, pcov = scipy.optimize.curve_fit(Gaussian, xydata, zdata, p0=guess )
		amp, xc, yc,  c = popt
		print("solution\n",amp, xc, yc, c)
		xCenter, yCenter = (popt[1], popt[2])
		radius = 2.35482*c / 2
		#yradius = 2.35*c
		#print("position and FWHM: ", xCenter, yCenter, xradius, yradius)
		print("position and FWHM: ", xCenter, yCenter, radius*2)
		if plot:
			matplotlib.pyplot.plot([xCenter, xCenter], [0, size], color="red")
			matplotlib.pyplot.plot([0, size], [yCenter, yCenter], color="red")
			# Draw a circle at the FWHM
			circle1 = matplotlib.pyplot.Circle((xCenter, yCenter), radius, color='r', fill=False)
			matplotlib.pyplot.gca().add_artist(circle1)
		
			# Plot the 3-d image
			matplotlib.pyplot.figure(plotGaussian.number)
			ax = plotGaussian.gca(projection="3d")
			surf = ax.plot_surface(x, y, z, linewidth=0, antialiased=True)
			matplotlib.pyplot.draw()

			# Plot a 2-D map of the solution
			matplotlib.pyplot.figure(fitFigure.number)
			fit = Gaussian(xydata, amp, xc, yc, c)
			fit.shape = (size, size)
			ax = matplotlib.pyplot.gca()
			im = ax.imshow(fit)
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.plot([xCenter, xCenter], [0, size], color="red")
			matplotlib.pyplot.plot([0, size], [yCenter, yCenter], color="red")
			circle1 = matplotlib.pyplot.Circle((xCenter, yCenter), radius, color='r', fill=False)
			matplotlib.pyplot.gca().add_artist(circle1)
			
			matplotlib.pyplot.show(block=True)

		apertureScaler = radius * 1.6
		apertures = CircularAperture(positions, r=apertureScaler)
		skyApertures = CircularAnnulus(positions, r_in=apertureScaler+10, r_out=apertureScaler+18)
		apers = [apertures, skyApertures]
		apertures.plot(color='red', lw=1.5, alpha=0.7)
		skyApertures.plot(color='blue', lw=1.5, alpha=0.7)
		error = numpy.sqrt(imageData)
		phot_table = aperture_photometry(imageData, apers, error=error)
		for col in phot_table.colnames:
			phot_table[col].info.format = '%.8g'  # for consistent table outputphot_table['aperture_sum'].info.format = '%.8g'
		bkg_mean = phot_table['aperture_sum_1']  / skyApertures.area
		bkg_sum = bkg_mean * apertures.area
		final_sum = phot_table['aperture_sum_0'] - bkg_sum
		phot_table['final_sum_err'] = numpy.sqrt(phot_table['aperture_sum_err_0']**2 + (phot_table['aperture_sum_err_1']/skyApertures.area)**2)
		phot_table['residual_aperture_sum'] = final_sum
		effective_gain = frame['exptime']  # seconds
		print(phot_table)
		for index, phot in enumerate(phot_table):
			measurement = {'JD': frame['JD'], 'ADUs': phot['residual_aperture_sum']}
			measurement['exptime'] = frame['exptime']
			measurement['countspersecond'] = phot['residual_aperture_sum'] / frame['exptime']
			measurement['error'] = phot['final_sum_err']
			measurement['extractposition'] = (phot['xcenter'].value, phot['ycenter'].value)
			targetList[index].addMeasurement(measurement)

		# Draw the image
		amplifiedImage = generallib.percentiles(imageData, 5, 95)
		matplotlib.pyplot.imshow(amplifiedImage)
		matplotlib.pyplot.gca().invert_yaxis()
		if frame['frame']==0:
			for i, p in enumerate(positions):
				matplotlib.pyplot.text(p[0]+10, p[1]+10, "%d"%i, color="blue", size="xx-large")
			matplotlib.pyplot.savefig("finder.png")
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(arg.pause)

		matplotlib.pyplot.clf()
		print("Frame number: {:d}".format(frame['frame']))
		if frame['frame']>=stopFrame-1: break


	photometryPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	


	colours = ['red', 'green', 'blue', 'purple', 'gray']
	for index, t in enumerate(targetList):
		measurements = t.getMeasurements()
		startDate = int(measurements[0]['JD'])
		xValues = [m['JD'] - startDate for m in measurements]
		yValues = [m['countspersecond'] for m in measurements]
		yErrors = [0 for m in measurements]
		print(yErrors)
		matplotlib.pyplot.errorbar(xValues, yValues, yerr=yErrors, color=colours[index], fmt = '.', ecolor=colours[index], capsize=0, label=index)
	matplotlib.pyplot.legend()
	matplotlib.pyplot.xlabel("JD - {:d}".format(startDate))
	matplotlib.pyplot.ylabel("counts per second")
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	
	differentialPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	
	comparison = targetList[1].getMeasurements()
	target = targetList[0].getMeasurements()
	startDate = int(target[0]['JD'])
	yValues = [t['countspersecond']/c['countspersecond'] for t, c, in zip(target, comparison)]
	xValues = [t['JD'] - startDate for t in target]
	matplotlib.pyplot.scatter(xValues, yValues, label="Target/comparison", s=3)
	matplotlib.pyplot.legend()
	matplotlib.pyplot.xlabel("JD - {:d}".format(startDate))
	matplotlib.pyplot.ylabel("Relative ADUs")
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	
	
	# write photometry to CSV file
	sample = targetList[0].getMeasurements()
	JDs = [s['JD'] for s in sample]
	exptimes = [s['exptime'] for s in sample]
	outputfile = open("photometry.csv", 'wt')
	outputfile.write("JD, exptime, ")
	for index in range(len(targetList)):
		outputfile.write("aperture_index, counts_per_second, extract_x, extract_y, ")
	outputfile.write("\n")

	for index in range(len(JDs)):
		outputfile.write("{:.6f}, {:.2f}, ".format(JDs[index], exptimes[index]))
		sys.stdout.write("{:.6f}, {:.2f}, ".format(JDs[index], exptimes[index]))
		for j, t in enumerate(targetList):
			m = t.getMeasurements()
			sys.stdout.write("{:d}, {:.2f}, ".format(j, m[index]['countspersecond']))
			sys.stdout.write("{:.2f}, {:.2f}, ".format(m[index]['extractposition'][0], m[index]['extractposition'][1]))
			outputfile.write("{:d}, {:.2f}, ".format(j, m[index]['countspersecond']))
			outputfile.write("{:.2f}, {:.2f}, ".format(m[index]['extractposition'][0], m[index]['extractposition'][1]))
			
		sys.stdout.write("\n")
		outputfile.write("\n")

	outputfile.close()
	input("Press enter to continue")
	
	sys.exit()
