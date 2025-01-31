import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from typing import Tuple, List, Dict
import utils
import pickle
import at
from copy import deepcopy

class OrbitCorrectionNN(nn.Module):
    def __init__(self, n_bpm_inputs: int, n_correctors: int, config_encoding_size: int):
        """
        Neural Network for orbit correction
        Args:
            n_bpm_inputs: Number of BPM readings (x,y for each BPM)
            n_correctors: Number of correctors (total of horizontal and vertical)
            config_encoding_size: Size of the configuration encoding based on physical parameters
        """
        super(OrbitCorrectionNN, self).__init__()
        
        # Separate networks for processing BPM readings and configuration
        self.bpm_network = nn.Sequential(
            nn.Linear(n_bpm_inputs, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU()
        )
        
        # Configuration network processes actual physical parameters
        self.config_network = nn.Sequential(
            nn.Linear(config_encoding_size, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU()
        )
        
        # Combined network for final prediction
        self.combined_network = nn.Sequential(
            nn.Linear(128 + 32, 64),  # 128 from BPM network + 32 from config network
            nn.ReLU(),
            nn.Linear(64, n_correctors)
        )
        
    def forward(self, bpm_data: torch.Tensor, config_data: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through the network
        Args:
            bpm_data: Tensor of BPM readings
            config_data: Tensor of configuration parameters
        """
        bpm_features = self.bpm_network(bpm_data)
        config_features = self.config_network(config_data)
        combined = torch.cat([bpm_features, config_features], dim=1)
        return self.combined_network(combined)

class OrbitCorrector:
    def __init__(self, base_ring, device='cpu'):
        """
        Orbit correction using ML
        Args:
            base_ring: Base accelerator lattice
            device: Device to run computations on ('cpu' or 'cuda')
        """
        self.base_ring = base_ring
        self.device = device
        self.seed_rings = {}  # Dictionary to store different seed configurations
        
        # Get number of BPMs and correctors from base ring
        self.bpm_indices = base_ring.get_uint32_index(at.Monitor)
        self.n_bpms = len(self.bpm_indices)
        
        # Get corrector counts
        self.hcm = utils.getCorrectorStrengths(base_ring, 'x')
        self.vcm = utils.getCorrectorStrengths(base_ring, 'y')
        self.n_correctors = len(self.hcm) + len(self.vcm)
        
        # Calculate configuration encoding size based on physical parameters we want to track
        self.config_encoding_size = self._calculate_config_encoding_size(base_ring)
        
        # Initialize network
        self.model = OrbitCorrectionNN(
            n_bpm_inputs=self.n_bpms * 2,  # x and y readings for each BPM
            n_correctors=self.n_correctors,
            config_encoding_size=self.config_encoding_size
        ).to(device)
        
        self.optimizer = optim.Adam(self.model.parameters())
        self.criterion = nn.MSELoss()
        
        # Store normalization parameters for configuration encoding
        self.config_norm_params = {}
    
    def _calculate_config_encoding_size(self, ring):
        """Calculate the size of configuration encoding based on physical parameters"""
        # Count the number of parameters we'll encode
        encoding_size = 0
        
        # Quadrupole strengths
        quad_indices = ring.get_uint32_index(at.Quadrupole)
        encoding_size += len(quad_indices)  # K1 values
        
        # RF parameters (if present)
        rf_indices = ring.get_uint32_index(at.RFCavity)
        if len(rf_indices) > 0:
            encoding_size += 2  # Voltage and Frequency

        # Dipole parameters
        dipole_indices = ring.get_uint32_index(at.Dipole)
        encoding_size += len(dipole_indices)  # Bend angles
        
        return encoding_size
    
    def _encode_configuration(self, ring) -> np.ndarray:
        """
        Extract physical parameters from ring configuration
        Returns normalized array of physical parameters
        """
        params = []
        
        # Extract Quadrupole strengths (K1 values)
        quad_indices = ring.get_uint32_index(at.Quadrupole)
        for idx in quad_indices:
            params.append(ring[idx].K)
            
        # Extract RF parameters if present
        rf_indices = ring.get_uint32_index(at.RFCavity)
        if len(rf_indices) > 0:
            params.extend([
                ring[rf_indices[0]].Voltage,
                ring[rf_indices[0]].Frequency
            ])
            
        # Extract Dipole parameters
        dipole_indices = ring.get_uint32_index(at.Dipole)
        for idx in dipole_indices:
            params.append(ring[idx].BendingAngle)
        
        return np.array(params)
    
    def _normalize_config_encoding(self, encoding: np.ndarray, is_training: bool = False) -> np.ndarray:
        """Normalize configuration encoding using mean and std"""
        if is_training:
            # During training, compute and store normalization parameters
            self.config_norm_params['mean'] = np.mean(encoding)
            self.config_norm_params['std'] = np.std(encoding)
        
        # Normalize using stored parameters
        return (encoding - self.config_norm_params['mean']) / (self.config_norm_params['std'] + 1e-8)
    
    def load_seed_configurations(self, seed_range=(1, 100)):
        """Load all seed configurations into memory"""
        for seed_num in range(seed_range[0], seed_range[1] + 1):
            seed_file = f'./lattices/seed{seed_num:03d}_pyAT_postRFcorrection'
            with open(seed_file, 'rb') as fid:
                self.seed_rings[seed_num] = pickle.load(fid)
    
    def prepare_input(self, bpm_readings: List[np.ndarray], ring) -> Tuple[torch.Tensor, torch.Tensor]:
        """Convert BPM readings and ring configuration to model input format"""
        # Process BPM readings
        x = np.concatenate(bpm_readings)
        bpm_tensor = torch.FloatTensor(x).to(self.device)
        
        # Process configuration
        config_encoding = self._encode_configuration(ring)
        normalized_encoding = self._normalize_config_encoding(config_encoding)
        config_tensor = torch.FloatTensor(normalized_encoding).to(self.device)
        
        return bpm_tensor, config_tensor
    
    def train(self, train_data, val_data, epochs=100, batch_size=32):
        """
        Train the neural network
        Args:
            train_data: List of (bpm_readings, ring, target_corrector_values) tuples
            val_data: Validation data in same format as train_data
            epochs: Number of training epochs
            batch_size: Batch size for training
        """
        train_losses = []
        val_losses = []
        
        # Compute normalization parameters from training data
        all_configs = np.vstack([self._encode_configuration(b[1]) for b in train_data])
        self._normalize_config_encoding(all_configs, is_training=True)
        
        for epoch in range(epochs):
            self.model.train()
            epoch_loss = 0
            batch_count = 0
            
            # Training loop
            for i in range(0, len(train_data), batch_size):
                batch = train_data[i:i + batch_size]
                bpm_inputs = torch.stack([self.prepare_input(b[0], b[1])[0] for b in batch])
                config_inputs = torch.stack([self.prepare_input(b[0], b[1])[1] for b in batch])
                
                self.optimizer.zero_grad()
                predicted_corrections = self.model(bpm_inputs, config_inputs)
                
                # Calculate loss using RMS difference
                batch_loss = torch.tensor(0.0, device=self.device, requires_grad=True)
                for j, (bpm_readings, ring, _) in enumerate(batch):
                    # Apply predicted corrections to ring
                    current_ring = deepcopy(ring)
                    n_hcm = len(self.hcm)
                    hcm_values = predicted_corrections[j, :n_hcm].detach().cpu().numpy()
                    vcm_values = predicted_corrections[j, n_hcm:].detach().cpu().numpy()
                    
                    current_ring = utils.setCorrectorStrengths(current_ring, 'x', hcm_values)
                    current_ring = utils.setCorrectorStrengths(current_ring, 'y', vcm_values)
                    
                    # Get new BPM readings
                    new_bpm_readings, new_trajectory = utils.getBPMreading(current_ring)
                    
                    # Calculate RMS difference loss
                    loss_value = torch.tensor(
                        sum(abs(utils.rms(new_bpm_readings) - utils.rms(new_trajectory))).astype('float32'),
                        device=self.device
                    )
                    batch_loss = loss_value + batch_loss
                
                batch_loss = batch_loss / len(batch)
                batch_loss.backward()
                self.optimizer.step()
                
                epoch_loss += batch_loss.item()
                batch_count += 1
            
            avg_train_loss = epoch_loss / batch_count
            train_losses.append(avg_train_loss)
            
            # Validation
            self.model.eval()
            val_loss = 0
            val_count = 0
            with torch.no_grad():
                for bpm_readings, ring, _ in val_data:
                    bpm_inputs, config_inputs = self.prepare_input(bpm_readings, ring)
                    bpm_inputs = bpm_inputs.unsqueeze(0)
                    config_inputs = config_inputs.unsqueeze(0)
                    predicted_corrections = self.model(bpm_inputs, config_inputs)
                    
                    # Apply corrections and calculate loss
                    current_ring = deepcopy(ring)
                    n_hcm = len(self.hcm)
                    hcm_values = predicted_corrections[0, :n_hcm].cpu().numpy()
                    vcm_values = predicted_corrections[0, n_hcm:].cpu().numpy()
                    
                    current_ring = utils.setCorrectorStrengths(current_ring, 'x', hcm_values)
                    current_ring = utils.setCorrectorStrengths(current_ring, 'y', vcm_values)
                    
                    new_bpm_readings, new_trajectory = utils.getBPMreading(current_ring)
                    val_loss += sum(abs(utils.rms(new_bpm_readings) - utils.rms(new_trajectory)))
                    val_count += 1
            
            avg_val_loss = val_loss / val_count
            val_losses.append(avg_val_loss)
            
            if (epoch + 1) % 10 == 0:
                print(f'Epoch {epoch+1}/{epochs}:')
                print(f'Training Loss: {avg_train_loss:.6f}')
                print(f'Validation Loss: {avg_val_loss:.6f}')
    
    def predict_corrections(self, bpm_readings: List[np.ndarray], ring) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict corrector values for given BPM readings and ring configuration
        Returns:
            Tuple of (horizontal_corrector_values, vertical_corrector_values)
        """
        self.model.eval()
        with torch.no_grad():
            bpm_inputs, config_inputs = self.prepare_input(bpm_readings, ring)
            bpm_inputs = bpm_inputs.unsqueeze(0)
            config_inputs = config_inputs.unsqueeze(0)
            predictions = self.model(bpm_inputs, config_inputs).cpu().numpy()[0]
            
            # Split predictions into horizontal and vertical corrections
            n_hcm = len(self.hcm)
            hcm_values = predictions[:n_hcm]
            vcm_values = predictions[n_hcm:]
            
            return hcm_values, vcm_values
    
    def apply_corrections(self, bpm_readings: List[np.ndarray], seed_num: int) -> None:
        """Apply ML-predicted corrections to the ring"""
        if seed_num not in self.seed_rings:
            raise ValueError(f"Seed configuration {seed_num} not loaded")
            
        ring = self.seed_rings[seed_num]
        hcm_values, vcm_values = self.predict_corrections(bpm_readings, ring)
        
        # Apply corrections to the specific seed configuration
        self.seed_rings[seed_num] = utils.setCorrectorStrengths(ring, 'x', hcm_values)
        self.seed_rings[seed_num] = utils.setCorrectorStrengths(self.seed_rings[seed_num], 'y', vcm_values) 