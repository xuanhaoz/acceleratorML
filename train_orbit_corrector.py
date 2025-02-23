import numpy as np
import torch
import pickle
from ml_orbit_correction import OrbitCorrector
import utils
from tqdm import tqdm
import logging
import os

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
    
    for seed_num in tqdm(range(seed_range[0], seed_range[1] + 1), 
                        desc="Generating training data", 
                        unit="seed"):
        # Load pre-correction configuration
        pre_seed_file = f'./matlab/seeds/seed{seed_num:d}_preCorrection_pyAT'
        post_seed_file = f'./matlab/seeds/seed{seed_num:d}_postCorrection_pyAT'
        
        try:
            with open(pre_seed_file, 'rb') as fid:
                pre_ring = pickle.load(fid)
            with open(post_seed_file, 'rb') as fid:
                post_ring = pickle.load(fid)
                
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
            
            # Store the data
            training_data.append((true_trajectory, initial_correctors, target_corrections))
            
        except FileNotFoundError:
            logger.warning(f"Skipping seed {seed_num} - files not found")
            continue
        except Exception as e:
            logger.error(f"Skipping seed {seed_num} - {e}")
            continue
    
    # Save to cache
    training_data = np.array(training_data, dtype=object)
    logger.info(f"Saving data cache to {cache_file}")
    np.savez(cache_file, training_data=training_data)

    return training_data





def main():
    # Load initial lattice as base configuration
    lattice_file = './matlab/seeds/seed1_preCorrection_pyAT'
    with open(lattice_file, 'rb') as fid:
        base_ring = pickle.load(fid)
    
    # Initialize orbit corrector
    logger.info("Initializing orbit corrector...")
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    corrector = OrbitCorrector(base_ring, device=device)
    
    # Generate training and validation data
    logger.info("Generating training data...")

    train_seeds = range(1, 801)
    val_seeds = range(801, 811)
    test_seeds = range(901, 1001)
    
    train_data = generate_training_data(seed_range=(train_seeds.start, train_seeds.stop - 1))
    
    # Log dataset sizes
    logger.info(f"Training samples: {len(train_data)}")
    logger.info(f"Validation samples: {len(val_seeds)}")
    
    # Train the model
    logger.info(f"Training model on {device}...")
    corrector.train(
        train_data=train_data,
        val_seeds=val_seeds,
        epochs=500,
        batch_size=128
    )

    logger.info("Testing model...")
    test_results = corrector.validate(test_seeds)


if __name__ == "__main__":
    main()