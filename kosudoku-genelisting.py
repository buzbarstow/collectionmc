#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# kosudoku-genelisting.py
# Created by Buz Barstow 2018-04-13
# Last modified by Buz Barstow 2018-11-09
#
# Code to make listing of gene loci in the Gluconobacter genome.
# ----------------------------------------------------------------------------------------------- #

from kosudoku.utils import ensure_dir, FindATandTAPositions, ImportGenBankSequence
from kosudoku.input import get_input
from kosudoku.feature import ImportFeatureArrayFromGenBank
from scipy import unique, setdiff1d, intersect1d
import re
import sys



# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'essentialGeneFileName', 'referenceGenomeGenBankFileName', \
'cdsOutputFileName', 'geneOutputFileName', 'genesWithoutCDSsFileName', 'cdssWithoutGenesFileName']

# argv = sys.argv

argv=['','inputdata/kosudoku-genelisting.inp']

inputParameterValues = get_input(argv, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
referenceGenomeGenBankFileName = inputParameterValues['referenceGenomeGenBankFileName']
essentialGeneFileName = inputParameterValues['essentialGeneFileName']

cdsOutputFileName = inputParameterValues['cdsOutputFileName']
geneOutputFileName = inputParameterValues['geneOutputFileName']
genesWithoutCDSsFileName = inputParameterValues['genesWithoutCDSsFileName']
cdssWithoutGenesFileName = inputParameterValues['cdssWithoutGenesFileName']
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Ensure that there are directories for output files
ensure_dir(cdsOutputFileName)
ensure_dir(geneOutputFileName)
ensure_dir(genesWithoutCDSsFileName)
ensure_dir(cdssWithoutGenesFileName)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)
sequence = ImportGenBankSequence(referenceGenomeGenBankFileName)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import the essential genes
essentialGeneHandle = open(essentialGeneFileName, 'r')
essentialGeneLines = essentialGeneHandle.readlines()
essentialGeneLociNames = []
for line in essentialGeneLines:
	essentialGeneLociNames.append(line.strip())
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
cdsFeatureArray = []
geneFeatureArray = []

for feature in featureArray:
	if feature.featureType == 'CDS':
		cdsFeatureArray.append(feature)
	elif feature.featureType == 'gene':
		geneFeatureArray.append(feature)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
cdsLocusAndFeatureNameArray = []
cdsLocusNameArray = []
cdsDispensableLocusNameArray = []

geneLocusAndFeatureNameArray = []
geneLocusNameArray = []
geneDispensableLocusNameArray = []


for feature in cdsFeatureArray:
	tagDict = feature.tagDict
	
	if 'locus_tag' in tagDict:
	
		locusTag = tagDict['locus_tag'][0]
		featureName = feature.featureName
		
		if locusTag in essentialGeneLociNames:
			essentiality = 'Essential'
		else:
			essentiality = 'Dispensable'
			cdsDispensableLocusNameArray.append(locusTag)
	
		cdsLocusAndFeatureNameArray.append([locusTag, featureName, essentiality])
		cdsLocusNameArray.append(locusTag)
		
for feature in geneFeatureArray:
	tagDict = feature.tagDict
	
	if 'locus_tag' in tagDict:
	
		locusTag = tagDict['locus_tag'][0]
		featureName = feature.featureName
		
		if locusTag in essentialGeneLociNames:
			essentiality = 'Essential'
		else:
			essentiality = 'Dispensable'
			geneDispensableLocusNameArray.append(locusTag)
	
		geneLocusAndFeatureNameArray.append([locusTag, featureName, essentiality])
		geneLocusNameArray.append(locusTag)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Figure out the CDSs that don't have genes, genes that don't have CDSs

genesWithoutCDSs = setdiff1d(geneLocusNameArray, cdsLocusNameArray)
cdssWithoutGenes = setdiff1d(cdsLocusNameArray, geneLocusNameArray)


uniqueGenes = unique(geneLocusNameArray)
uniqueDispensableGenes = unique(geneDispensableLocusNameArray)
uniqueCDSs = unique(cdsLocusNameArray)
uniqueDispensableCDSs = unique(cdsDispensableLocusNameArray)

# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
# Write out data

cdsFileHandle = open(cdsOutputFileName, 'w')

for line in cdsLocusAndFeatureNameArray:
	outputStr = str(line[0]) + ',' + str('"' + line[1] + '"') + ',' + str(line[2]) + '\n'
	cdsFileHandle.write(outputStr)

cdsFileHandle.close()



geneFileHandle = open(geneOutputFileName, 'w')

for line in geneLocusAndFeatureNameArray:
	outputStr = str(line[0]) + ',' + str(line[1]) + ',' + str(line[2]) + '\n'
	geneFileHandle.write(outputStr)

geneFileHandle.close()



genesWithoutCDSsHandle = open(genesWithoutCDSsFileName, 'w')

for line in genesWithoutCDSs:
	outputStr = line + '\n'
	genesWithoutCDSsHandle.write(outputStr)

genesWithoutCDSsHandle.close()



cdssWithoutGenesHandle = open(cdssWithoutGenesFileName, 'w')

for line in cdssWithoutGenes:
	outputStr = line + '\n'
	cdssWithoutGenesHandle.write(outputStr)

cdssWithoutGenesHandle.close()


# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
# Print some output data

print("Unique genes: " + str(len(uniqueGenes)))
print("Unique CDSs: " + str(len(uniqueCDSs)))
print("Unique dispensable genes: " + str(len(uniqueDispensableGenes)))
print("Unique dispensable CDSs: " + str(len(uniqueDispensableCDSs)))
print("Genes without CDSs: " + str(len(genesWithoutCDSs)))
print("CDSs without genes: " + str(len(cdssWithoutGenes)))
# ----------------------------------------------------------------------------------------------- #


