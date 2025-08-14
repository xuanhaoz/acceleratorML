# AcceleratorML Project Code Summary

## Quick Overview

This project focuses on developing machine learning solutions for orbit correction in particle accelerators. Both supervised and reinforcement learning methods have been explored, along with several deep learning architectures. While the focus remains the deployment of the current model on the synchrotron, we should explore other approaches more thoroughly if time allows. 

## Key Dependencies

- **pyAT** (Accelerator Toolbox): Core accelerator physics simulation
- **PyTorch**: Deep learning framework
- **TorchRL**: Reinforcement learning framework
- **MATLAB SC Toolkit**: Seed and reference correction generation

## Project Structure and Components

### 1. Core Python Implementation

#### Original Example Skeleton (`main.py`)
- Demonstrates basic orbit correction workflow and usage of key functions
- Loads lattice data from pickle files (seed001_pyAT_postRFcorrection)
- Shows how to:
  - Get BPM readings using `utils.getBPMreading()`
  - Extract corrector magnet strengths
  - Apply ML-generated corrections
  - Calculate loss metrics

#### Utility Functions (`utils.py`)
Core accelerator physics functions including:

**Key Functions:**
- `getBPMreading(ring, makePlot=0)`: Get BPM readings with measurement errors
  - Supports both closed orbit and single-shot tracking modes (don't worry about those)
  - Handles particle losses and injection errors
  - Includes BPM calibration errors and noise
- `getCorrectorStrengths(ring, plane)`: Extract horizontal/vertical corrector values
- `setCorrectorStrengths(ring, plane, setpoints)`: Apply corrector settings
- `checkClosedOrbitExists(ring)`: Verify orbit stability (internal use)

### 2. Machine Learning Components

#### Neural Network Orbit Corrector (`ml_orbit_correction.py`)

**OrbitCorrectionNN Class:**
- Deep neural network (4096→2048→1024→n_correctors)
- Uses BatchNorm, ReLU activations, and Xavier initialization
- Input: BPM trajectory data concatenated with initial corrector values
- Output: Corrector strength adjustments

**OrbitCorrector Class:**
- Complete ML pipeline for orbit correction
- Features data normalization using StandardScaler
- Training loop with validation on separate seed datasets
- Prediction with proper denormalization

**Key Methods:**
- `train()`: Main training loop with batching and validation
- `predict_corrections()`: Generate corrector adjustments for new trajectories
- `validate()`: Test model performance on validation seeds

#### Training Pipeline (`train_orbit_corrector.py`)
- Automated training data generation from pre/post correction seed pairs
- Multiprocessing support (16 workers) for parallel data processing
- Caches relevant data extracted from seeds to avoid regenerating on each run

### 3. Reinforcement Learning Implementation
May be relevant in the future if we want to retry this approach
#### RL Training Notebook (`rl_training.ipynb`)
Complete RL solution using PPO (Proximal Policy Optimization):

**Environment (`LatticeEnv`):**
- Parallelized simulation environment (64 parallel processes)
- Custom TorchRL environment for accelerator physics
- State: BPM readings (360 BPMs × 2 coordinates)
- Action: Corrector magnet adjustments (792 correctors)
- Reward: Negative RMS orbit error

**PPO Agent:**
- Actor-Critic architecture with separate networks
- Policy network outputs Gaussian distributions for continuous control
- Value network estimates state values for advantage calculation
- GAE (Generalized Advantage Estimation) for improved training stability

### 4. Advanced Data Processing

#### Diffusion Dataset (`diffusion_dataset.py`)
Experimental diffusion model approach:
- DDPM (Denoising Diffusion Probabilistic Model) scheduler
- Generates training chains by adding noise to corrector settings
- Multi-step denoising process for gradual orbit correction
- Preprocessing for large-scale diffusion model training

### 5. Data and Lattice Configuration

#### Lattice Files (`lattices/` directory)
- Not currently used
- 100+ pre-generated lattice configurations with different error seeds
- Post-RF correction lattices for each seed
- Variety of accelerator configurations (AS2v225, ASSR variants)

#### MATLAB Reference Implementation (`matlab/` directory)
Not necessary to know, as details have been abstracted away for us.
For those who are curious to know more though, perhaps ask Frank for better explanations.
**Key Components:**
- **SC Toolkit**: Complete accelerator simulation framework
- **runOrbitCorrection_AS2.m**: Reference orbit correction algorithm
- **Response Matrices**: Pre-calculated correction matrices
- **Lattice Definitions**: MATLAB lattice definitions for different configurations

**Features:**
- SVD-based orbit correction with regularization
- Progressive correction with multiple alpha values
- Support for both orbit and pseudo-orbit modes
- Sextupole ramping for commissioning scenarios

### 6. Machine Interface (`machine/` directory)

#### Real Machine Integration
- `example_beam_loss_monitor_program`: Interface to actual accelerator control system
- EPICS integration for real-time machine control
- Beam loss monitoring capabilities
- Production-ready code for deployment

#### Process Variable Lists (`PVlist`)
- Comprehensive list of machine process variables
- Maps simulation parameters to real machine controls


## Current Status and Results

### Working Components
- **Data pipeline**: Automated generation and caching  
- **ML training**: Neural network successfully trained and evaluated
- **RL implementation**: PPO agent trained and validated  
- **Baseline comparison**: MATLAB reference implementation functional  

### Experimental Components
- **Diffusion models**: Alternative generative approach (Leo: this is a crackpot idea that I had one day and never finished implementing. Ask me if you want to know more)
- **Advanced architectures**: Exploring transformer and CNN-based models


## Future Directions

### Immediate Development
- **Real machine deployment**: Deploy trained models on actual accelerator
- **Performance benchmarking**: Comprehensive comparison with traditional methods
- **Robustness testing**: Validation under various machine conditions

### Research Extensions
- **Alternative architectures**: Transformers/CNNs are likely less prone to overfitting
- **Surpassing SVD**: RL will likely be required to surpass performance of traditional methods
- **Online learning**: Continuous model improvement during operation

## Technical Notes

### Data Formats
- **Lattice files**: Python pickle format (pyAT rings)
- **Training data**: NumPy compressed format (.npz)
- **Models**: PyTorch state dictionaries

### Computational Requirements
- **Training**: GPU very much recommended
- **Inference**: CPU sufficient for real-time operation with current model
