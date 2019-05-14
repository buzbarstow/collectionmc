import xml.etree.ElementTree as ET
import pdb

importedCoords = {}

tree = ET.parse('out.xml')
root = tree.getroot()
importedCoordsList = root.findall('coord')

for coord in importedCoordsList:
		
	coordinate = int(coord.attrib['coord'])
	loci = coord.findall('locus')
# 	pdb.set_trace()
	
	importedCoordsKeys = importedCoords.keys()
	
	if coordinate not in importedCoordsKeys:
		importedCoords[coordinate] = []
	
	for locus in loci:
		locusName = locus.attrib['locus']
		essentiality = locus.attrib['essentiality']
		importedCoords[coordinate].append([locusName, essentiality])
