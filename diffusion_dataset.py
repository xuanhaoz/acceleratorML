# %%
import torch
import at
import os
import io
from tqdm import tqdm
import utils
import numpy as np
import math
from sklearn.preprocessing import StandardScaler
import scipy.io as sio



# %% [markdown]
# so what we want to do is save the dataset in a format that allows for more flexibility during training later. for each seed, we will save a few chains of noising steps. each chain will save: the bpm readings at each noising step and the cm readings at each noising step. each will be a separate tensor, so e.g. indexing i+1 will yield the prev step readings. should save each chain as a different filename inside the same npz.

# %%
"""
class DDPM_Scheduler(nn.Module):
    def __init__(self, num_time_steps: int=1000):
        super().__init__()
        self.beta = torch.linspace(1e-4, 0.02, num_time_steps, requires_grad=False)
        alpha = 1 - self.beta
        self.alpha = torch.cumprod(alpha, dim=0).requires_grad_(False)

    def forward(self, t):
        return self.beta[t], self.alpha[t]
"""
    
def make_ddpm_scheduler(time_steps=1000, min_beta=1e-4, max_beta=0.02):
    global out_scheduler
    def out_scheduler(t):
        beta = np.linspace(min_beta, max_beta, num=time_steps)
        alpha = 1 - beta
        return beta[t], alpha[t]
    return out_scheduler


def check_file(file):
    try:
        matfile = sio.loadmat(file)
    except:
        return False

    post_cor = matfile.get('postCorrection')
    if len(post_cor) != 4154:
        return False
    else:
        return True
    

# %%
def worker_process(args):
    save_name, nd_save_name, noise_scheduler, scaler, path_str, max_steps, chains_per_seed = args
    
    id = int(''.join(filter(str.isdigit, path_str)))

    try:
        no_correction = at.load_mat(path_str, use="preCorrection")
        init_ring = at.load_mat(path_str, use="postCorrection")

        # saving the pre correction bpm readings just in case. no need to save
        # cm as they are all initialised with 0
        pre_bpm, _ = utils.getBPMreading(no_correction)
        bpm, _ = utils.getBPMreading(no_correction)

        hcm = utils.getCorrectorStrengths(init_ring, 'x')
        vcm = utils.getCorrectorStrengths(init_ring, 'y') 
        cm = np.append(hcm, vcm)
        np.savez(nd_save_name.format(x=id), pre_bpm=pre_bpm, post_bpm=bpm, cm=cm) 
    except:
        # error with seed, skipping...
        return

    # creating chains
    for i in range(chains_per_seed):

        # making one chain of corrections
        ring = init_ring.deepcopy()
        cms = []
        bpms = []
        prev_ring = init_ring.deepcopy()
        t = 0
        while t < max_steps:
            try:
                # get current bpm readings
                bpm, _ = utils.getBPMreading(ring)
            except:
                # unable to get bpms, revert to previous step
                ring = prev_ring.deepcopy()
                ring = add_noise(ring, hcm, vcm, t - 1, noise_scheduler, scaler)
                continue
            bpms.append(bpm)

            # get current corrector magnet settings
            hcm = utils.getCorrectorStrengths(ring, 'x')
            vcm = utils.getCorrectorStrengths(ring, 'y') 
            cms.append(np.append(hcm, vcm))

            prev_ring = ring.deepcopy()
            # add noise
            ring = add_noise(ring, hcm, vcm, t, noise_scheduler, scaler)
            t += 1

        cms = np.stack(cms, 0)
        bpms = np.stack(bpms, 0)
        np.savez(save_name.format(x=id, y=i), bpms=bpms, cms=cms)


def add_noise(ring, hcm, vcm, t, noise_scheduler, scaler: StandardScaler):
    cm = np.append(hcm, vcm)
    hcm_sep = len(hcm)
    noise = np.random.normal(size=cm.shape)
    noise = noise * scaler.scale_
    beta, alpha = noise_scheduler(t)
    noised_cm = math.sqrt(alpha) * cm + math.sqrt(beta) * noise

    hcm = noised_cm[:hcm_sep]
    vcm = noised_cm[hcm_sep:]

    ring = utils.setCorrectorStrengths(ring, 'x', hcm)
    ring = utils.setCorrectorStrengths(ring, 'y', vcm)

    return ring


# %%
from tqdm.contrib.concurrent import process_map
import tqdm
class DiffusionPreProcessor():
    def __init__(self, num_workers: int, noise_cheduler, 
                 save_name: str, nd_save_name: str, max_steps=50, chains_per_seed=16):
        self.num_workers = num_workers
        self.noise_scheduler = noise_cheduler
        self.save_name = save_name
        self.max_steps = max_steps
        self.chains_per_seed = chains_per_seed
        self.nd_save_name = nd_save_name
        
    def process_seeds(self, seed_dir, scaler, id_range=None):
        if id_range is None:
            id_range = (0, len(os.listdir(seed_dir)))
        seed_files = os.listdir(seed_dir)[id_range[0]:id_range[1]]
        worker_args = [(self.save_name, self.nd_save_name, self.noise_scheduler, 
                        scaler, seed_dir + seed_file, self.max_steps, self.chains_per_seed) for seed_file in seed_files]
        # tqdm_class=tqdm.notebook.tqdm
        process_map(worker_process, worker_args, max_workers=self.num_workers)


# %%
noise_scheduler = make_ddpm_scheduler(min_beta=0.000000000001, max_beta=0.1, time_steps=128) #0.002

if __name__ == "__main__":
    seed_dir = "./matlab/seeds/"
    num_workers = 4
    save_name = "./diffusion_data/seed{x}_chain{y}.spz"
    max_steps = 128
    chains_per_seed = 2
    nd_save_name = "./single_step_data/seed{x}.npz"

    cached_data = np.load("data_cache_0_15999.npz", allow_pickle=True)
    corrected = cached_data['training_data'][:, 2]
    corrected2 = np.stack(corrected)

    scaler = StandardScaler(with_mean = False)
    scaler.fit(corrected2)

    data_preprocessor = DiffusionPreProcessor(num_workers, noise_scheduler, save_name, nd_save_name, max_steps, chains_per_seed)

    # %%
    data_preprocessor.process_seeds(seed_dir, scaler)


