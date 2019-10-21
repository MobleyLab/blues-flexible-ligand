# This script is for finding total number of moves accepted as the simulatio progresses.
# Output: #movesAccepted #movesProposed
# Usage: python3 findAcceptanceIteration.py inputFilename outputFilename 

import sys

# opening output file
output = open( sys.argv[2], 'w' )


# opening input file
with open(sys.argv[1]) as fid:
   line = fid.readline()
   while line:
    prevLine = line
    line = fid.readline()
    if "Performing " in line:
       splitted = prevLine.split()
       output.write( splitted[3]+"\t"+str(int(splitted[4])+1)+"\n" )

output.close()
