#!/sw/bin/python3.6

# ----------------------------------------------------------------------------------------------- #
# kosudoku-findhimar.py
# Created by Buz Barstow 2018-04-01
# Last modified by Buz Barstow 2018-10-23
#
# Code to find percentage of coding sequences in genome that can have the himar AT/TA
# insertion site
# ----------------------------------------------------------------------------------------------- #

from kosudoku.utils import ensure_dir, FindATandTAPositions, ImportGenBankSequence
from kosudoku.feature import ImportFeatureArrayFromGenBank
from kosudoku.input import get_input


import re
import sys
import pdb


# ----------------------------------------------------------------------------------------------- #
inputParameters = [\
'referenceGenomeGenBankFileName']

# argv = sys.argv

argv=['','inputdata/kosudoku-findhimar.inp']

inputParameterValues = get_input(argv, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
referenceGenomeGenBankFileName = inputParameterValues['referenceGenomeGenBankFileName']
# ----------------------------------------------------------------------------------------------- #




# ----------------------------------------------------------------------------------------------- #
# Parse the input data
featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)
sequence = ImportGenBankSequence(referenceGenomeGenBankFileName)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Find the features that are coding regions and the AT/TA positions

atRegex = re.compile('(at|ta)', re.IGNORECASE)

cdsFeatureArray = []
cdsATorTANumberArray = []

for feature in featureArray:
	if feature.featureType == 'CDS':
		startIndex = feature.startCoord-1
		endIndex = feature.endCoord-1
		
		featureSequence = sequence[startIndex:endIndex+1]
		feature.sequence = featureSequence
		
		featureATandTAPositions = []
		i = 0
		
		
		while i < len(featureSequence) - 1:
			atMatch = atRegex.match(sequence[i:i+2])
		
			if atMatch != None:
				featureATandTAPositions.append(i+1)
		
			i += 1
		
		
		feature.ATandTAPositions = featureATandTAPositions
		
		cdsATorTANumberArray.append(len(featureATandTAPositions))
		cdsFeatureArray.append(feature)
		
	

numberFeaturesWithATorTA = 0
totalFeatures = 0
for feature in cdsFeatureArray:
	totalFeatures += 1
	if len(feature.ATandTAPositions) > 0:
		numberFeaturesWithATorTA += 1

percentageFeaturesWithATorTA = (float(numberFeaturesWithATorTA)/float(totalFeatures))*100.
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Find the total number of AT and TA positions in the genome
		
ATandTAPositions = []
i = 0
while i < len(sequence) - 1:
	atMatch = atRegex.match(sequence[i:i+2])
	
	if atMatch != None:
		ATandTAPositions.append(i+1)
	
	i += 1

totalNumberOfATandTAPositions = len(ATandTAPositions)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
print("Total number of AT and TA positions: " + str(totalNumberOfATandTAPositions))
print("Total number of features: " + str(totalFeatures))
print("Total number of features with AT or TA: " + str(numberFeaturesWithATorTA))
print("Percentage of features with AT or TA: " + str(percentageFeaturesWithATorTA))
# ----------------------------------------------------------------------------------------------- #

