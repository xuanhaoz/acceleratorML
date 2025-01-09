# utility functions
# created 11/Dec/2024, F. Zhang
# major changes:
# optimized performance and added multiprocessing utility functions, L. Qiao

import at
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from scipy.stats import truncnorm

import time
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
def getBPMreading(ring,makePlot=0,use_guess=False):
    """ get BPM reading from tracking and add measurement errors
        Leo - This has been modified to only call find_orbit6() once, which
        should roughly halve the run time. The use_guess param specifies 
        whether we should accelerate convergence by using the last orbit 
        as a guess.

    Args:
        ring (lattice array): lattice array
        makePlot (int, optional): flag for plotting BPM readings. Defaults to 0.

    Returns:
        B (len(bpm),2): stacked array for BPM reading in x and y with errors
        T (len(ring),2): stacked array for relative position between beam trajectory and element centre
    """
    start = time.time()

    if not hasattr(getBPMreading, "guess") or use_guess == False:
        orbit, T = at.find_orbit6(ring,range(len(ring)))
    else:
        orbit, T = at.find_orbit6(ring, range(len(ring)), guess=getBPMreading.guess)
    getBPMreading.guess = orbit

    closedOrbitExists =  ~np.any(np.isnan(orbit))
    check_closed = time.time()

    bpm = ring.get_uint32_index(at.Monitor)

    inds = time.time()

    if ~closedOrbitExists:
        # track in single shot mode
        #
        print("single shot mode")
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
        #T = at.find_orbit6(ring,range(len(ring)))
        #T = T[1]
        pass
    
    find_orbit = time.time()

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

    find_offset = time.time()

    # calculate effects on BPM reading due to BPM errors
    # NOT INCLUDED: TbT noise
    #
    BPMcalError = np.array([ele.CalError for ele in ring[bpm]])
    BPMnoise    = np.array([ele.NoiseCO for ele in ring[bpm]])

    get_noise = time.time()

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

    print(f"""getBPMreading: 
          checking closed took {check_closed - start} secs
          getting BPM indices took {inds - check_closed} secs
          finding orbit took {find_orbit - inds} secs
          finding offset took {find_offset - find_orbit} secs
          getting noise took {get_noise - find_offset} secs""")
    return B,T

# ------------------------------------------------
#
def getCorrectorStrengths(ring,plane):
    # plane: x = horizontal, y = vertical
    #
    chv = ring.get_uint32_index(at.Corrector)
    if plane == 'x':
        # focusing sextupoles - hori. correcotr (HCM)
        #
        sf = ring.get_bool_index("SF*")
        sf = np.where(sf)[0]
        ords = sorted(np.append(chv,sf))
    elif plane == 'y':
        # defocusing sextupoels - ver. corr. (VCM)
        #
        sd = ring.get_bool_index("SD*")
        sd = np.where(sd)[0]
        ords = sorted(np.append(chv,sd))

    monitor = np.zeros(len(ords))
    for n,ele in enumerate(ring[ords]):
        if plane == 'x':
            polynom = ele.PolynomB[0]
        elif plane == 'y':
            polynom = ele.PolynomA[0]

        monitor[n] = polynom

    return monitor

# ------------------------------------------------
#
def setCorrectorStrengths(ring,plane,setpoints,in_place=False):
    """Leo - I've added an option to forgo the deepcopy operation and simply
       modify the ring in place. This option is significantly faster, so pass
       in_place=True whenever possible.
    """
    # plane: 0 = horizontal, 1 = vertical
    #
    start = time.time()
    chv = ring.get_uint32_index(at.Corrector)
    if plane == 'x':
        # focusing sextupoles - hori. correcotr (HCM)
        #
        sf = ring.get_bool_index("SF*")
        sf = np.where(sf)[0]
        ords = sorted(np.append(chv,sf))
    elif plane == 'y':
        # defocusing sextupoels - ver. corr. (VCM)
        #
        sd = ring.get_bool_index("SD*")
        sd = np.where(sd)[0]
        ords = sorted(np.append(chv,sd))
    
    inds = time.time()

    assert len(ords) == len(setpoints), f"Ords have lenght {len(ords)} but setpoints have length {len(setpoints)}"
    
    if not in_place:
        newRing = deepcopy(ring)
    else:
        newRing = ring
        
    copied = time.time()

    for n,ele in enumerate(newRing[ords]):
        setpoint = setpoints[n]
        if plane == 'x':
            ele.PolynomB[0] = setpoint
        elif plane == 'y':
            ele.PolynomA[0] = setpoint

    set_elems = time.time()

    #print(f"""setCorrectorStrengths: 
    #      getting indices took {inds - start} secs
    #      deepcopying took {copied - inds} secs
    #      setting elements took {set_elems - copied} secs""")
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


# ------------------------------------------------
# Multiprocessing helper functions
import torch

NAN_COST = 1
def worker_set_corrector(cm, ring, hcm_sep, device):
    hcm = cm[:hcm_sep]
    vcm = cm[hcm_sep:]
    ring = setCorrectorStrengths(ring, 'x', hcm, True)
    ring = setCorrectorStrengths(ring, 'y', vcm, True)
    
    bpm, traj = getBPMreading(ring)
    
    cost = torch.tensor(traj, dtype=torch.float32, device=device)
    cost = torch.linalg.vector_norm(cost, ord=2, dim=-1)

    cost = torch.nan_to_num(cost, NAN_COST)
    cost = torch.sum(torch.abs(cost))
    bpm = torch.nan_to_num(torch.tensor(bpm, dtype=torch.float32, device=device), NAN_COST)
    return ring, cost, bpm


def update_ring(ring, cm, hcm_sep, nan_cost):
    hcm = cm[:hcm_sep]
    vcm = cm[hcm_sep:]
    ring = setCorrectorStrengths(ring, 'x', hcm, True)
    ring = setCorrectorStrengths(ring, 'y', vcm, True)

    bpm, traj = getBPMreading(ring)

    # currently using the raw euclidean distance to element centres as loss,
    # maybe should square?
    cost = np.linalg.norm(traj, axis=-1)
    cost = np.nan_to_num(cost, copy=False, nan=nan_cost)
    cost = np.sum(np.abs(cost))

    bpm = np.nan_to_num(bpm, copy=False, nan=nan_cost)
    return cost, bpm

from multiprocessing.connection import Connection
    
def worker_process(ring: at.lattice.lattice_object.Lattice, conn: Connection):
    hcm = getCorrectorStrengths(ring, 'x')
    vcm = getCorrectorStrengths(ring, 'y') 
    hcm_sep = len(hcm)
    cm = np.append(hcm, vcm)
    init_cm = cm.copy()
    NAN_COST = 1

    while True:
        # we receive a tuple of (command: str, data) from the parent process
        command, data = conn.recv()

        if command == "step":
            cm += data
            out = update_ring(ring, cm, hcm_sep, NAN_COST)
            conn.send(out)
        
        elif command == "reset":
            cm = init_cm.copy()
            out = update_ring(ring, cm, hcm_sep, NAN_COST)
            conn.send(out)
            
        elif command == "ping":
            conn.send("received")