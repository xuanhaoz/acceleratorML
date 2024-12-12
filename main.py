# main function
# created 11/Dec/2024, F. Zhang
# major changes:
#

import at
import pickle
import numpy as np
from copy import deepcopy
import utils

latticeFile = './lattices/seed001_pyAT_postRFcorrection'
with open(latticeFile,'rb') as fid:
    ring0 = pickle.load(fid)

# get current BPM readings
# returns:
#   B0: [x,y] reading at every BPM
#   T0: [x,y] relative offset between beam trajectory and element centre at every element
#
[B0, T0] = utils.getBPMreading(ring0)

# get current corrector settings
#
HCM = utils.getCorrectorStrengths(ring0,'x')
VCM = utils.getCorrectorStrengths(ring0,'y')
CM = np.append(HCM,VCM)

# inputs - BPM readings
#
flattenB = np.concatenate(B0)

# ML training...
# output should be of dimension:
#
outputShape = len(CM)
MLoutput = np.random.randn(outputShape) * 1e-5

# update CM values and get new BPM readings
#
HCM1 = MLoutput[:len(HCM)]
VCM1 = MLoutput[len(HCM):]

newRing = deepcopy(ring0)
newRing = utils.setCorrectorStrengths(newRing,'x',HCM1)
newRing = utils.setCorrectorStrengths(newRing,'y',VCM1)

[B1,T1] = utils.getBPMreading(newRing)

# potential loss value:
#
loss = sum(abs(utils.rms(B1) - utils.rms(T1)))
print(f"Loss value: {loss}")



