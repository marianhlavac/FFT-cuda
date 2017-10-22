import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 3:
    print('Wrong count of arguments.')
    sys.exit()

inputdatafile = sys.argv[1]
outputdatafile = sys.argv[2]

inputchart = plt.plotfile(inputdatafile)
plt.title('Input data')
plt.xlabel('time')
plt.ylabel('value')

outputchart = plt.plotfile(outputdatafile, (0, 1), delimiter=',')
plt.title('Output data (FFT)')

plt.show()