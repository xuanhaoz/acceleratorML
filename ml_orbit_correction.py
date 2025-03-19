import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from typing import List
from tqdm import tqdm
import pickle
import at

import utils
from sklearn.preprocessing import StandardScaler

class OrbitCorrectionNN(nn.Module):
    def __init__(self, n_elements: int, n_correctors: int, dropout_rate: float = 0.5):
        """Neural Network for orbit correction using layers"""
        super(OrbitCorrectionNN, self).__init__()

        # Calculate input size including initial corrector values
        input_size = n_elements * 2 + n_correctors

        # Define network architecture
        self.network = nn.Sequential(
            # First layer
            nn.Linear(input_size, 4096),
            nn.BatchNorm1d(4096),
            nn.ReLU(),

            # Second layer
            nn.Linear(4096, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),

            # Third layer
            nn.Linear(2048, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU(),

            # Output layer
            nn.Linear(1024, n_correctors),
        )

        self._initialize_weights()

    def _initialize_weights(self):
        """Initialize network weights using Xavier initialization"""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    def forward(self, x):
        return self.network(x)


class OrbitCorrector:
    def __init__(self, base_ring, device='cpu', dropout_rate=0.5, weight_decay=1e-5):
        """Initialize orbit corrector with ML model"""
        self.base_ring = base_ring
        self.device = device
        self.weight_decay = weight_decay

        # Initialize scalers
        self._init_scalers()

        # Setup model parameters
        self._setup_model_params()

        # Initialize model and optimizer
        self.model = OrbitCorrectionNN(
            n_elements=self.n_true_trajectory_inputs,
            n_correctors=self.n_correctors,
            dropout_rate=dropout_rate
        ).to(device)

        self.optimizer = optim.Adam(self.model.parameters(), weight_decay=weight_decay)
        self.criterion = nn.MSELoss()

    def _init_scalers(self):
        """Initialize data scalers"""
        self.trajectory_scaler = StandardScaler()
        self.corrector_scaler = StandardScaler()
        self.initial_corrector_scaler = StandardScaler()

    def _setup_model_params(self):
        """Setup model parameters from base ring"""
        bpm_readings, true_trajectory = utils.getBPMreading(self.base_ring)
        self.n_true_trajectory_inputs = len(bpm_readings)

        self.hcm = utils.getCorrectorStrengths(self.base_ring, 'x')
        self.vcm = utils.getCorrectorStrengths(self.base_ring, 'y')
        self.n_correctors = len(self.hcm) + len(self.vcm)

    def fit_scalers(self, train_data):
        """Fit the scalers on training data"""
        trajectories = np.vstack([
            data[0].reshape(-1, 2) if data[0].ndim == 1 else data[0]
            for data in train_data
        ])
        corrections = np.vstack([data[2] for data in train_data])
        initial_corrections = np.vstack([data[1] for data in train_data])

        self.trajectory_scaler.fit(trajectories)
        self.corrector_scaler.fit(corrections)
        self.initial_corrector_scaler.fit(initial_corrections)

    def prepare_input(self, true_trajectory: np.ndarray, initial_correctors: np.ndarray) -> torch.Tensor:
        """
        Prepare trajectory and initial corrector data for model input
        Args:
            true_trajectory: Array of shape (n_elements, 2) with x,y positions
            initial_correctors: Array of initial corrector values
        """
        # Reshape trajectory to 2D array if needed
        if true_trajectory.ndim == 1:
            true_trajectory = true_trajectory.reshape(-1, 2)

        # Transform using fitted scalers
        normalized_traj = self.trajectory_scaler.transform(true_trajectory)
        normalized_init_corr = self.initial_corrector_scaler.transform(
            initial_correctors.reshape(1, -1)
        )

        # Concatenate and flatten for model input
        x = np.concatenate([normalized_traj.flatten(), normalized_init_corr.flatten()])

        return torch.FloatTensor(x).to(self.device)

    def prepare_target(self, target_corrections: np.ndarray) -> torch.Tensor:
        """Normalize target corrections using fitted scaler"""
        normalized_corrections = self.corrector_scaler.transform(
            target_corrections.reshape(1, -1) if target_corrections.ndim == 1
            else target_corrections
        )
        return torch.FloatTensor(normalized_corrections).to(self.device)


    def train(self, train_data, val_seeds, epochs=100, batch_size=32):
        """Train the neural network with normalized data"""
        # First fit the scalers on training data
        self.fit_scalers(train_data)

        train_losses = []
        val_losses = []

        for epoch in range(epochs):
            self.model.train()
            epoch_loss = 0
            batch_count = 0

            # Training loop
            for i in range(0, len(train_data), batch_size):
                batch = train_data[i:i + batch_size]

                # Prepare normalized inputs and targets
                inputs = torch.stack([
                    self.prepare_input(b[0], b[1]) for b in batch
                ])
                target_corrections = torch.stack([
                    self.prepare_target(b[2])[0] for b in batch
                ])

                self.optimizer.zero_grad()
                predicted_corrections = self.model(inputs)

                loss = self.criterion(predicted_corrections, target_corrections)
                loss.backward()
                self.optimizer.step()

                epoch_loss += loss.item()
                batch_count += 1

            avg_train_loss = epoch_loss / batch_count
            train_losses.append(avg_train_loss)

            # Validation
            if(epoch+1)%50==0:
                self.model.eval()
                with torch.no_grad():
                    self.validate(val_seeds)

            if (epoch + 1) % 10 == 0:
                print(f'Epoch {epoch+1}/{epochs}:')
                print(f'Training Loss: {avg_train_loss:.6f}')

        return train_losses, val_losses

    def validate(self, val_seeds):
        """Test trained model on test seeds and calculate losses"""
        results = []
        pbar = tqdm(val_seeds, desc="Testing model")
        for seed_num in pbar:
            # Load pre and post correction rings
            lattice_file = f"./matlab/seeds/seed{seed_num:d}.mat"

            try:
                pre_ring = at.load_mat(lattice_file, check=False, use="preCorrection")
                post_ring = at.load_mat(lattice_file, check=False, use="postCorrection")

                # Get initial state
                [B0, T0] = utils.getBPMreading(pre_ring)
                initial_rms = utils.rms(np.concatenate(T0))
                initial_loss = utils.rms(B0)

                # Get expected metrics
                target_hcm = utils.getCorrectorStrengths(post_ring, 'x')
                target_vcm = utils.getCorrectorStrengths(post_ring, 'y')
                pre_ring = utils.setCorrectorStrengths(pre_ring, 'x',target_hcm)
                pre_ring = utils.setCorrectorStrengths(pre_ring, 'y',target_vcm)
                [B1, T1] = utils.getBPMreading(pre_ring)
                expected_rms = utils.rms(np.concatenate(T1))
                expected_loss = utils.rms(B1)


                # Get current trajectory
                initial_hcm = utils.getCorrectorStrengths(pre_ring, 'x')
                initial_vcm = utils.getCorrectorStrengths(pre_ring, 'y')
                initial_correctors = np.concatenate([initial_hcm, initial_vcm])

                # Get model predictions for corrector settings
                predicted_corrections = self.predict_corrections(T0, initial_correctors)

                # Apply predicted corrections
                pre_ring = utils.setCorrectorStrengths(pre_ring, 'x',
                                                           predicted_corrections[:len(self.hcm)])
                pre_ring = utils.setCorrectorStrengths(pre_ring, 'y',
                                                           predicted_corrections[len(self.hcm):])

                # Measure new state
                [B_new, T_new] = utils.getBPMreading(pre_ring)

                new_rms = utils.rms(np.concatenate(T_new))
                new_loss = utils.rms(B_new)


                results.append({
                    'seed': seed_num,
                    'rms_improvement': ((initial_rms - new_rms) / initial_rms) * 100,
                    'expected rms_improvement': ((initial_rms - expected_rms) / initial_rms) * 100,
                    'loss_improvement': ((initial_loss - new_loss) / initial_loss) * 100,
                    'expected loss_improvement': ((initial_loss - expected_loss) / initial_loss) * 100,
                })

            except Exception as e:
                print(f"Error processing seed {seed_num}: {str(e)}")
                continue

        print("total rms improvement:" + str(np.mean([r['rms_improvement'] for r in results])) + "%")
        print("expected rms improvement:" + str(np.mean([r['expected rms_improvement'] for r in results])) + "%")
        print("total loss improvement:" + str(np.mean([r['loss_improvement'] for r in results])) + "%")
        print("expected loss improvement:" + str(np.mean([r['expected loss_improvement'] for r in results])) + "%")

        return results

    def predict_corrections(self, true_trajectory: List[np.ndarray], initial_correctors: np.ndarray) -> np.ndarray:
        """Predict corrector values with proper denormalization"""
        self.model.eval()
        with torch.no_grad():
            true_trajectory_inputs = self.prepare_input(true_trajectory, initial_correctors)
            true_trajectory_inputs = true_trajectory_inputs.unsqueeze(0)
            predictions = self.model(true_trajectory_inputs)

            # Denormalize predictions
            denormalized_predictions = self.corrector_scaler.inverse_transform(
                predictions.cpu().numpy()
            )
            return denormalized_predictions[0]
