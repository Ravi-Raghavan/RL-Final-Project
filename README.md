# Satellite Attitude Takeover Control Simulation

This repository contains MATLAB simulation code for optimizing control strategies for satellite attitude takeover, based on the linear quadratic differential game approach outlined in the paper "Linear quadratic differential game approach for attitude takeover control of failed spacecraft" by Yuan Chai, Jianjun Luo, Nan Han, and Jianfeng Xie.

## Overview

The simulation implements a Reinforcement Learning algorithm to control the attitude of a failed spacecraft by utilizing microsatellites. It leverages Nash Differential Games to derive optimal control strategies and applies Lyapunov iterations for convergence of the control solutions.

## Prerequisites

- MATLAB
- Control System Toolbox for MATLAB (for functions like `are` and `lyap`)

## Files Description

- `Reinforcement_Learning_LQ_Nash_Strategies.m`: The main script that sets up the simulation environment, including the inertia matrix, control gains, and initial conditions, and then runs the simulation using the Ordinary Differential Equation (ODE) solver.

## Usage

1. Clone the repository to your local machine.
2. Open MATLAB and navigate to the cloned repository directory.
3. Run the `Reinforcement_Learning_LQ_Nash_Strategies.m` script to start the simulation.

## Methodology

The code first defines the inertia matrix of a hypothetical failed spacecraft and computes its eigenvalues and eigenvectors. It then derives the control matrices for the microsatellites involved in the takeover. The simulation sets up the linear quadratic regulator (LQR) problem and iteratively solves the coupled algebraic Riccati equations (CARE) using the Lyapunov iterations method.

## Visualization

The code plots the following results from the simulation:
- Attitude angles over time (`γ`, `θ`, `ψ`)
- Angular velocity components over time (`w_x`, `w_y`, `w_z`)
- Control torques provided by each microsatellite over time (`u_1`, `u_2`, `u_3`)

## Contributors

- **Yuan Chai** - Original paper author and theoretical framework provider.
- **Jianjun Luo** - Original paper author and algorithm designer.
- **Nan Han** - Original paper author and provided mathematical proofs.
- **Jianfeng Xie** - Original paper author and simulation advisor.

## Citation

Please cite the following paper when using this code for academic purposes:

Chai, Y., Luo, J., Han, N., & Xie, J. (Year). Linear quadratic differential game approach for attitude takeover control of failed spacecraft. Chai, Y., Luo, J., Han, N., & Xie, J. (2020). Linear quadratic differential game approach for attitude takeover control of failed spacecraft. Acta Astronautica, 175, 142–154. https://doi.org/10.1016/j.actaastro.2020.04.023


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details.

## Acknowledgments

This project owes its existence to the hard work and dedication of its contributors. We extend our deepest gratitude to the following individuals who have played a pivotal role in the development and success of this simulation:

- The authors of the referenced paper for their groundbreaking work on the subject.
- **Ravi Raghavan** - For his valuable contributions to the core algorithms and simulation design.
- **Rahul Hegde** - For his expertise in system modeling and insightful analysis throughout the development process.
- **Prayag Patel** - For his innovative solutions in reinforcement learning application and optimization strategies.
- **James Chan** - For his meticulous approach to code integrity, testing, and validation efforts.

Their collective efforts have significantly advanced the project's goals and have set a high standard for excellence in collaborative development.


## Contact

For any queries regarding the simulation code, please raise an issue in the repository or contact the maintainers.
