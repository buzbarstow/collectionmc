#!/sw/bin/python3.6

# ----------------------------------------------------------------------------------------------- #
# kosudoku-collectionmc
# Created by Buz Barstow 2016-05-06
# Last modified by Buz Barstow 2018-11-09

# Monte Carlo code that simulates picking of a random progenitor library and determines gene
# knockout coverage. 
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
from kosudoku.montecarlo import SimulateMultiplePickings, SimulatePicking, PoissonEstimateOfGenesHit
from kosudoku.output import generateOutputMatrixWithHeaders, writeOutputMatrix
from kosudoku.input import get_input
from kosudoku.utils import ensure_dir

import sys
from numpy import ones, float
from matplotlib.pyplot import ion, ioff, show, figure, plot, grid, xlabel, ylabel, legend
# ------------------------------------------------------------------------------------------------ #


# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file
# argv = sys.argv

argv = ['','inputdata/kosudoku-collectionmc-ReutrophaH16Nottingham-gene.inp']

inputParameters = ['numberOfTrials', 'maxMutants', 'outputFileName', \
'transposonCoordToFeatureIndexFileName']

inputParameterValues = get_input(argv, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
numberOfTrials = int(inputParameterValues['numberOfTrials'])
maxMutants = int(inputParameterValues['maxMutants'])
transposonCoordToFeatureIndexFileName = \
inputParameterValues['transposonCoordToFeatureIndexFileName']
outputFileName = inputParameterValues['outputFileName']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
ensure_dir(outputFileName)
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Do the simulation for the maximum non essential set

[iAxis, averageFeatureHitCount, sdFeatureHitCount, featureHitCountUpperBound, \
featureHitCountLowerBound, noUniqHittableFeatures] = \
SimulateMultiplePickings(transposonCoordToFeatureIndexFileName, numberOfTrials, maxMutants)

noUniqHittableFeaturesArray = ones(len(iAxis), float)*noUniqHittableFeatures
noUniqHittableFeaturesArray99 = ones(len(iAxis), float)*noUniqHittableFeatures*0.99
noUniqHittableFeaturesArray95 = ones(len(iAxis), float)*noUniqHittableFeatures*0.95

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Do a poisson simulation as well
poissonUniqueGenesHit = PoissonEstimateOfGenesHit(iAxis, noUniqHittableFeatures)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
headers = ['Mutants', 'AvUniqFeatures', 'SDUniq', 'FeaturesUpper', 'FeaturesLower']
vectorList = [iAxis, averageFeatureHitCount, sdFeatureHitCount, \
featureHitCountUpperBound, featureHitCountLowerBound]
oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
writeOutputMatrix(outputFileName, oMatrix)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output simulation data
figure()
plot(iAxis, averageFeatureHitCount, c='purple', label='Monte Carlo estimate')
plot(iAxis, featureHitCountUpperBound, c='blue')
plot(iAxis, featureHitCountLowerBound, c='blue')

plot(iAxis, noUniqHittableFeaturesArray, c='black')
plot(iAxis, noUniqHittableFeaturesArray99, c='black')
plot(iAxis, noUniqHittableFeaturesArray95, c='black')

plot(iAxis, poissonUniqueGenesHit, c='red', label='Poisson estimate')
grid()
xlabel('Mutants Picked')
ylabel('Non-essential Genes Hit')
legend(loc=4)
show()
# ------------------------------------------------------------------------------------------------ #
