import numpy as np
import torch
from copy import deepcopy
import at
import pickle
from ml_orbit_correction import OrbitCorrector
import utils
from tqdm import tqdm

def apply_static_errors(ring):
    """Apply static errors as described in the paper"""
    # Element-to-element misalignment: 15 µm rms
    element_misalignment = np.random.normal(0, 15e-6, len(ring))
    
    # Girder center misalignment: 15 µm rms
    girder_misalignment = np.random.normal(0, 15e-6, len(ring))
    
    # Magnets rotation: 150 µrad rms
    magnet_rotation = np.random.normal(0, 150e-6, len(ring))
    
    # Apply misalignments to ring elements
    ring_with_errors = deepcopy(ring)
    
    # Apply the errors to the ring elements
    for i, element in enumerate(ring_with_errors):
        if hasattr(element, 'Offset'):
            current_offset = np.array(element.Offset)
            element.Offset = current_offset + np.array([element_misalignment[i], girder_misalignment[i]])
        if hasattr(element, 'Roll'):
            element.Roll = element.Roll + magnet_rotation[i]
    
    return ring_with_errors

def apply_dynamic_errors(ring):
    """Apply dynamic errors as described in the paper"""
    # Quadrupoles' Δk1/k1: ±1%
    quad_indices = ring.get_uint32_index(at.Quadrupole)
    ring_with_errors = deepcopy(ring)
    
    for idx in quad_indices:
        variation = np.random.uniform(-0.01, 0.01)  # ±1%
        ring_with_errors[idx].K = ring_with_errors[idx].K * (1 + variation)
    
    return ring_with_errors

def generate_training_data(seed_range=(1, 100), samples_per_seed=10):
    """
    Generate training data for the ML model using multiple seed configurations
    Args:
        seed_range: Tuple of (start_seed, end_seed) inclusive
        samples_per_seed: Number of samples to generate for each seed configuration
    """
    training_data = []
    
    for seed_num in tqdm(range(seed_range[0], seed_range[1] + 1), 
                        desc="Generating training data", 
                        unit="seed"):
        # Load seed configuration
        seed_file = f'./lattices/seed{seed_num:03d}_pyAT_postRFcorrection'
        with open(seed_file, 'rb') as fid:
            seed_ring = pickle.load(fid)
        
        # Apply static errors once per seed
        ring_with_static = apply_static_errors(seed_ring)
        
        for _ in range(samples_per_seed):
            # Apply dynamic errors
            current_ring = apply_dynamic_errors(deepcopy(ring_with_static))
            
            # Get initial BPM readings
            bpm_readings, _ = utils.getBPMreading(current_ring)
            
            # Store the data (we don't need target values anymore)
            training_data.append((bpm_readings, current_ring, None))
    
    return training_data

def evaluate_model(corrector, test_seeds):
    """
    Evaluate the model on multiple seed configurations
    Args:
        corrector: Trained OrbitCorrector instance
        test_seeds: List of seed numbers to test on
    """
    results = {}
    
    for seed_num in test_seeds:
        # Get the seed configuration
        if seed_num not in corrector.seed_rings:
            print(f"Skipping seed {seed_num} - not loaded")
            continue
            
        test_ring = deepcopy(corrector.seed_rings[seed_num])
        
        # Apply errors
        test_ring = apply_static_errors(test_ring)
        test_ring = apply_dynamic_errors(test_ring)
        
        # Get initial readings
        initial_readings, _ = utils.getBPMreading(test_ring)
        initial_rms = utils.rms(np.concatenate(initial_readings))
        
        # Apply ML-based correction
        corrector.apply_corrections(initial_readings, seed_num)
        
        # Get final readings
        final_readings, _ = utils.getBPMreading(corrector.seed_rings[seed_num])
        final_rms = utils.rms(np.concatenate(final_readings))
        
        # Store results
        results[seed_num] = {
            'initial_rms': initial_rms,
            'final_rms': final_rms,
            'improvement': (initial_rms - final_rms) / initial_rms * 100,
            'configuration': test_ring  # Store the actual configuration for analysis
        }
    
    return results

def main():
    # Load initial lattice as base configuration
    lattice_file = './lattices/seed001_pyAT_postRFcorrection'
    with open(lattice_file, 'rb') as fid:
        base_ring = pickle.load(fid)
    
    # Initialize orbit corrector
    print("Initializing orbit corrector...")
    device = 'cuda' if torch.cuda.is_available() else 'mps'
    corrector = OrbitCorrector(base_ring, device=device)
    
    # Load all seed configurations
    print("Loading seed configurations...")
    corrector.load_seed_configurations()
    
    # Generate training and validation data
    print("Generating training data...")
    train_seeds = range(1, 81)  # Use seeds 1-80 for training
    val_seeds = range(81, 91)   # Use seeds 81-90 for validation
    test_seeds = range(91, 101) # Use seeds 91-100 for testing
    
    train_data = generate_training_data(seed_range=(train_seeds.start, train_seeds.stop - 1), 
                                      samples_per_seed=50)
    val_data = generate_training_data(seed_range=(val_seeds.start, val_seeds.stop - 1), 
                                    samples_per_seed=10)
    
    # Train the model
    print(f"Training model on {device}...")
    corrector.train(
        train_data=train_data,
        val_data=val_data,
        epochs=100,
        batch_size=32
    )
    
    # Save the trained model
    torch.save({
        'model_state_dict': corrector.model.state_dict(),
        'optimizer_state_dict': corrector.optimizer.state_dict(),
    }, 'orbit_corrector_model.pth')
    print("Model saved to orbit_corrector_model.pth")
    
    # Evaluate the model
    print("\nEvaluating model on test seeds...")
    results = evaluate_model(corrector, test_seeds)
    
    # Print results
    print("\nTest Results:")
    for seed_num, metrics in results.items():
        print(f"\nSeed {seed_num}:")
        print(f"Initial RMS orbit deviation: {metrics['initial_rms']:.6f}")
        print(f"Final RMS orbit deviation: {metrics['final_rms']:.6f}")
        print(f"Improvement: {metrics['improvement']:.2f}%")

if __name__ == "__main__":
    main() 