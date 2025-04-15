import numpy as np
import torch
import pickle
from ml_orbit_correction import OrbitCorrector
import utils
from tqdm import tqdm
import logging
import os
import at
import multiprocessing
import glob

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def generate_training_data(seed_range=(1, 100), cache_dir='./data_cache'):
    """
    Generate training data using pre and post correction seeds with caching
    Args:
        seed_range: Tuple of (start_seed, end_seed) inclusive
        cache_dir: Directory to store cached data
    """
    # Create cache directory if it doesn't exist
    os.makedirs(cache_dir, exist_ok=True)
    
    # Generate cache filename based on seed range
    cache_file = os.path.join(cache_dir, f'data_cache_{seed_range[0]}_{seed_range[1]}.npz')
    
    # Try to load from cache first
    if os.path.exists(cache_file):
        logger.info(f"Loading cached data from {cache_file}")
        cached_data = np.load(cache_file, allow_pickle=True)
        return cached_data['training_data']
    
    logger.info(f"Generating new training data for seeds {seed_range[0]} to {seed_range[1]}")
    training_data = []
    
    # Load pre-correction configuration
    with multiprocessing.Pool(16) as p:
        training_data = list(tqdm(p.imap(pool_worker, [(f'./matlab/seeds/seed{seed_num:d}.mat', seed_num, logger)
                                                        for seed_num in range(seed_range[0], seed_range[1] + 1)]),
                                    total=seed_range[1] + 1 - seed_range[0],
                                    desc="Generating training data",
                                    unit="seed"))

    training_data = [i for i in training_data if not isinstance(i[0], type(None))]
    # Save to cache
    training_data = np.array(training_data, dtype=object)
    logger.info(f"Saving data cache to {cache_file}")
    np.savez(cache_file, training_data=training_data)

    return training_data

def pool_worker(in_tuple):
    seed_file, seed_num, logger = in_tuple
    try:
        pre_ring = at.load_mat(seed_file, check=False, use="preCorrection")
        post_ring = at.load_mat(seed_file, check=False, use="postCorrection")
            
        # Get BPM readings from pre-correction ring
        bpm_readings, true_trajectory = utils.getBPMreading(pre_ring)
        # Get target corrector values from pre-correction ring
        initial_hcm = utils.getCorrectorStrengths(pre_ring, 'x')
        initial_vcm = utils.getCorrectorStrengths(pre_ring, 'y')
        initial_correctors = np.concatenate([initial_hcm, initial_vcm])
        
        # Get target corrector values from post-correction ring
        target_hcm = utils.getCorrectorStrengths(post_ring, 'x')
        target_vcm = utils.getCorrectorStrengths(post_ring, 'y')
        target_corrections = np.concatenate([target_hcm, target_vcm])
        return bpm_readings, initial_correctors, target_corrections
    except FileNotFoundError:
        logger.warning(f"Skipping seed {seed_num} - files not found")
        return None, None, None
    except Exception as e:
        logger.error(f"Skipping seed {seed_num} - {e}")
        return None, None, None

def main():
    # Load initial lattice as base configuration
    lattice_file = "./matlab/seeds/seed0.mat"
    base_ring = at.load_mat(lattice_file, check=False, use="preCorrection")
    
    # Initialize orbit corrector
    logger.info("Initializing orbit corrector...")
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    corrector = OrbitCorrector(base_ring, device=device)
    
    # Generate training and validation data
    logger.info("Generating training data...")

    train_seeds = range(16000)
    val_seeds = range(16000, 16100)
    test_seeds = range(18000, 18100)
    
    train_data = generate_training_data(seed_range=(train_seeds.start, train_seeds.stop - 1))
    
    # Log dataset sizes
    logger.info(f"Initial training samples: {len(train_data)}")
    logger.info(f"Validation samples: {len(val_seeds)}")
    
    # Train the model with progressive augmentation
    logger.info(f"Training model on {device} with progressive augmentation...")
    corrector.train(
        train_data=train_data,
        val_seeds=val_seeds,
        epochs=100000,
        batch_size=1024,
        augment_paitience=2000
    )

    logger.info("Testing model...")
    # Load the best model before final testing
    logger.info("Loading best model for final testing...")
    best_model_files = glob.glob('saved_models/best_loss_model_*.pt')
    best_model_file = sorted(best_model_files, 
                           key=lambda x: float(x.split('improvement_')[1].split('pct.pt')[0]), 
                           reverse=True)[0]
    
    corrector.load_model_and_scalers(best_model_file)
    improvement = float(best_model_file.split('improvement_')[1].split('pct.pt')[0])
    logger.info(f"Loaded best model: {best_model_file} with {improvement:.2f}% improvement")
    
    test_results = corrector.validate(test_seeds)


if __name__ == "__main__":
    main()