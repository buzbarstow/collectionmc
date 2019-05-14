import xml.etree.ElementTree as ET

# ------------------------------------------------------------------------------------------------ #
def ExportListForXML(listToExport, delimeter=','):
	
	outputStr = ''
	
	i = 0
	while i < len(listToExport):
		outputStr += str(listToExport[i])
		if i < len(listToExport) - 1:
			outputStr += delimeter
		i += 1
	
	return outputStr
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
# ------------------------------------------------------------------------------------------------ #


coords = {\
1:[['locus1', 'dispensable'], ['locus2', 'dispensable']], \
50:[['locus2', 'dispensable'], ['locus4', 'dispensable']], \
100:[['locus5', 'dispensable'], ['locus6', 'dispensable']], \
201:[['locus10', 'essential']]}


coordListRoot = ET.Element('coordList')

coordKeys = coords.keys()

for coordKey in coordKeys:
	coordSubElement = ET.SubElement(coordListRoot, 'coord')
	coordSubElement.set('coord', str(coordKey))
	
	locusArray = coords[coordKey]
	
	for locus in locusArray:
		locusSubElement = ET.SubElement(coordSubElement, 'locus')
		locusSubElement.set('locus', locus[0])
		locusSubElement.set('essentiality', locus[1])
	


indent(coordListRoot)
coordListTree = ET.ElementTree(coordListRoot)
coordListTree.write('out.xml')





