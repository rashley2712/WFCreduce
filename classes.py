import json, numpy

class apertureClass:
	def __init__(self, position, peak):
		self.rootPosition = position
		self.peak = peak
		self.tracker = False
		self.target = False

	def __str__(self):
		return "peak {:.0f} position ({:.2f}, {:.2f})".format(self.peak, self.rootPosition[0], self.rootPosition[1])

	def toJSON(self):
		jsonObject = {}
		jsonObject['position'] = self.rootPosition
		jsonObject['peak'] = self.peak
		jsonObject['tracker'] = self.tracker
		jsonObject['target'] = self.target
		return jsonObject

	def setAsTarget(self, number):
		self.target= number

		
class apertureDB:
	def __init__(self):
		self.allApertures = []	

	def add(self, aperture):
		self.allApertures.append(aperture)

	def sort(self):
		def getPeak(e):
			return e.peak
		self.allApertures.sort(key = getPeak, reverse=True)

	def enableTrackers(self, n=10):
		for i, a in enumerate(self.allApertures):
			a.tracker = True
			if i>n: break

	def save(self, filename="apertures.json"):
		outputfile = open(filename, 'wt')
		jsonObject = []
		for a in self.allApertures:
			jsonObject.append(a.toJSON())
		json.dump(jsonObject, outputfile, indent=4)
		outputfile.close()

	def load(self, filename="apertures.json"):
		inputfile = open(filename, "rt")
		jsonObject = json.load(inputfile)
		for o in jsonObject:
			print(o)
			aperture = apertureClass(o['position'], o['peak'])
			aperture.tracker = o['tracker']
			aperture.target = o['target']
			self.allApertures.append(aperture)

	def makeCatalog(self):
		cat = []
		for a in self.allApertures:
			if a.tracker:
				cat.append([a.rootPosition[0], a.rootPosition[1]])
		return cat
	
	def getTargets(self):
		targetSubset = []
		for a in self.allApertures:
			if a.target: targetSubset.append(a)
		for t in targetSubset:
			print(t)
		targetSubset.sort(key=lambda x: x.target, reverse=False)
		for t in targetSubset:
			print(t)
		return targetSubset

	def getTargetByNumber(self, number):
		for aper in self.allApertures:
			if aper.target==number: return aper
		return None

	def match(self, x, y, number):
		radius = 5
		closest = 1E8
		bestMatch = -1
		for index, aper in enumerate(self.allApertures):
			distance = numpy.sqrt((aper.rootPosition[0]-x)**2 + (aper.rootPosition[1]-y)**2)
			if distance<closest:
				closest = distance
				bestMatch=index
		if closest>radius:
			print("No match within %d pixels."%radius)
			return None
		print("closest match is aperture: %s at %.1f pixels."%(self.allApertures[bestMatch], closest))
		self.allApertures[bestMatch].setAsTarget(number)
		return self.allApertures[bestMatch]

class frameDB:
	def __init__(self):
		self.allFrames = []

	def add(self, frame):
		self.allFrames.append(frame)

	def save(self, filename="frames.json"):
		outputfile = open(filename, 'wt')
		json.dump(self.allFrames, outputfile, indent=4)
		outputfile.close()

	def load(self, filename="frames.json"):
		inputfile = open(filename, "rt")
		self.allFrames = json.load(inputfile)
		print(self.allFrames)
		inputfile.close()

class target():
	def __init__(self, id):
		self.id = id
		self.name = "unknown"
		self.measurements = []

	def addMeasurement(self, measurement):
		self.measurements.append(measurement)

	def getMeasurements(self):
		return self.measurements