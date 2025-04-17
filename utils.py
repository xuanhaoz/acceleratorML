# utility functions
# created 11/Dec/2024, F. Zhang
# major changes:
#

import at
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from scipy.stats import truncnorm

# ------------------------------------------------
#
def rms(A):
    return np.sqrt(np.nanmean(A**2,axis=0))

# ------------------------------------------------
#
def checkClosedOrbitExists(ring):
    T = at.find_orbit6(ring)
    return ~np.any(np.isnan(T[0]))
    

# ------------------------------------------------
#
def getBPMreading(ring,makePlot=0):
    """ get BPM reading from tracking and add measurement errors

    Args:
        ring (lattice array): lattice array
        makePlot (int, optional): flag for plotting BPM readings. Defaults to 0.

    Returns:
        B (len(bpm),2): stacked array for BPM reading in x and y with errors
        T (len(ring),2): stacked array for relative position between beam trajectory and element centre
    """
    closedOrbitExists = checkClosedOrbitExists(ring)

    bpm = ring.get_uint32_index(at.Monitor)
    if ~closedOrbitExists:
        # track in single shot mode
        #
        nTurns = 1
        nShots = 10
        sigCutoff = 2
        injectionError = [1e-4, 1e-5, 1e-4, 1e-5, 1e-4, 1e-4] # [m, rad, m, rad, rel., m]

        Zin = np.zeros((6,nShots))
        for n,error in enumerate(injectionError):
            generator = truncnorm(-sigCutoff,sigCutoff,loc=0,scale=error)
            Zin[n,:] = generator.rvs(nShots)

        # Zout: [6 x Nshots x refpts x nTurns]
        #
        [Zout,trackParam,trackData] = ring.track(Zin,nTurns,refpts=range(len(ring)),losses=1)
        
        # at.track only changes x coord to nan if particle is lost
        # change all coords to nan for lost particles
        #
        IsNan = np.isnan(Zout)
        Zout[:,IsNan[0,:,:,0],:] = np.nan

        # average over all shots and transpose to get to same shape as closed orbit mode
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Do we want this nanmean here? This will bias bpm readings to shots that have small deviation
        # and ignore shots that have particle losses
        #
        T = np.nanmean(Zout[:,:,:,0],axis=1).T            
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        #
    else:
        # closed orbit mode
        #
        T = at.find_orbit6(ring,range(len(ring)))
        T = T[1]

    if np.isnan(T).any():
        raise Exception("Complete particle loss")

    # calculate relative different between beam trajectory and element centre
    # adding in element offset and roll errors
    #
    offset = np.zeros((len(ring),2))
    roll = np.zeros(len(ring))

    for n,ele in enumerate(ring):
        if hasattr(ele,'Offset'):
            offset[n,:] = np.array(ele.Offset) + np.array(ele.SupportOffset)

        if hasattr(ele,'Roll'):
            roll[n] = ele.Roll + ele.SupportRoll

    # only need x,y coords
    #
    T = T[:,[0,2]]

    # add offset and roll error
    #
    T[:,0] = np.cos(roll) * T[:,0] - np.sin(roll) * T[:,1]
    T[:,1] = np.sin(roll) * T[:,0] + np.cos(roll) * T[:,1]
    T = T - offset

    # calculate effects on BPM reading due to BPM errors
    # NOT INCLUDED: TbT noise
    #
    BPMcalError = np.array([ele.CalError for ele in ring[bpm]])
    BPMnoise    = np.array([ele.NoiseCO for ele in ring[bpm]])

    # add BPM calibration error and noise
    #
    B = T[bpm,:]
    B = B * (1 + BPMcalError)
    B = B + BPMnoise

    if makePlot:
        spos = at.get_s_pos(ring,range(len(ring)))
        fig = plt.figure()
        gs = fig.add_gridspec(2,hspace=0)
        axs = gs.subplots(sharex=True)

        axs[0].plot(spos[bpm],B[:,0]*1e6)
        axs[0].plot(spos,T[:,0]*1e6)
        axs[0].set_ylabel(f'x [$\mu$m]')
        axs[0].grid()

        axs[1].plot(spos[bpm],B[:,1]*1e6)
        axs[1].plot(spos,T[:,1]*1e6)
        axs[1].set_ylabel(f'y [$\mu$m]')
        axs[1].grid()

        plt.xlabel('s [m]')
        plt.legend(['BPMs','All eles'])
        plt.tight_layout()

    return B,T

# ------------------------------------------------
#
def getCorrectorStrengths(ring,plane,ver="assr"):
    # plane: x = horizontal, y = vertical
    #
    chv = ring.get_uint32_index(at.Corrector)
    match plane:
        case 'x':
            # focusing sextupoles - hori. correcotr (HCM)
            #
            sf = ring.get_bool_index("SF*")
            sf = np.where(sf)[0]
            ords = sorted(np.append(chv,sf))
            if ver == "assr":
                ords = sorted(sf)
        case 'y':
            # defocusing sextupoels - ver. corr. (VCM)
            #
            sd = ring.get_bool_index("SD*")
            sd = np.where(sd)[0]
            ords = sorted(np.append(chv,sd))
            if ver == "assr":
                ords = sorted(sd)

    monitor = np.zeros(len(ords))
    for n,ele in enumerate(ring[ords]):
        match plane:
            case 'x':
                polynom = ele.PolynomB[0]

            case 'y':
                polynom = ele.PolynomA[0]

        monitor[n] = polynom

    return monitor

# ------------------------------------------------
#
def setCorrectorStrengths(ring,plane,setpoints,ver="assr"):
    # plane: 0 = horizontal, 1 = vertical
    #
    chv = ring.get_uint32_index(at.Corrector)
    match plane:
        case 'x':
            # focusing sextupoles - hori. correcotr (HCM)
            #
            sf = ring.get_bool_index("SF*")
            sf = np.where(sf)[0]
            ords = sorted(np.append(chv,sf))
            if ver == "assr":
                ords = sorted(sf)
        case 'y':
            # defocusing sextupoels - ver. corr. (VCM)
            #
            sd = ring.get_bool_index("SD*")
            sd = np.where(sd)[0]
            ords = sorted(np.append(chv,sd))
            if ver == "assr":
                ords = sorted(sd)

    assert len(ords) == len(setpoints), f"Ords have lenght {len(ords)} but setpoints have length {len(setpoints)}"
    newRing = deepcopy(ring)

    for n,ele in enumerate(newRing[ords]):
        setpoint = setpoints[n]
        match plane:
            case 'x':
                ele.PolynomB[0] = setpoint
            case 'y':
                ele.PolynomA[0] = setpoint
        
        if ele.PassMethod == 'CorrectorPass':
            ele.KickAngle[0] = ele.PolynomB[0]
            ele.KickAngle[1] = ele.PolynomA[0]

    return newRing
        
# ------------------------------------------------
#
def initialiseNewSeed(ring,errorScale=1):
    # ring: initial pyAT ring lattice with errors already initialised 
    #   from Simulated Commissioning toolkit
    # errorScale: controls the scaling of random errors applied to the lattice
    #   0 -> inf. [float]
    #
    
    """
    newRing = deepcopy(ring)

    # define errors
    #
    errors = {}
    errors['roll'] = 100e-6
    errors['Dx'] = 30e-6
    errors['Dy'] = 30e-6

    for key,val in errors.items():
        errors[key] = val*errorScale

    """
    pass
