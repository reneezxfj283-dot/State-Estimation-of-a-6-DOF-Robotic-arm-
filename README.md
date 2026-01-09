# State-Estimation-of-a-6-DOF-Robotic-arm-
This repository provides the MATLAB simulation code accompanying a study on state estimation for a six-degree-of-freedom (6-DOF) robotic arm. The focus is on estimating joint positions and joint velocities under multiple models jumping due to payload mutation, using advanced variational Bayesian and multi-model estimation techniques.
The repository includes implementations of:
  The proposed method based on variational Bayesian (IEE)
  Interacting Multiple Model (IMM) filtering
  Classical Kalman filtering (KF)
Repository Structure:
├── rob_Main.m
│   └─ Main entry point for simulation and performance comparison
│
├── IEE_algorithm.m
│   └─ The proposed algorithm based on variational Bayesian 
│
├── IMM_algorithm.m
│   └─ Interacting Multiple Model filtering implementation
│
├── KF_algorithm.m
│   └─ Standard Kalman filter implementation
│
├── Dynamic.m
│   └─ Dynamic model of the robotic arm after linearization and discretization
│
├── linear.m
│   └─ Linearized system model
│
├── Rob_Model.m
│   └─ Setting the coefficients of the state-space equations of robotic arm
│
├── rob_Para.m
│   └─ Physical and model parameters of the robotic arm
│
├── DH.m
│   └─ Denavit–Hartenberg parameter definition
│
├── extractDynamicsComponents.m
│   └─ Extraction of inertia, Coriolis, and gravity-related terms
│
├── Gauss.m
│   └─ Gaussian noise generation utility
│
└── README.md
Requirements:
  MATLAB R2024b
  No external toolboxes are required
