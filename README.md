**Newton-Raphson Load Flow Analysis (MATLAB)**

This repository contains a suite of MATLAB functions and scripts designed to solve the Load Flow (Power Flow) Problem in electrical power systems using the Newton-Raphson (NR) Method.

Load flow analysis is fundamental to power system planning and operation, providing the steady-state voltages, angles, and power flows for a given set of generation and load conditions.

**🚀 Features**

Support for Multiple Bus Types: Handles Slack, PV (Generator), and PQ (Load) buses.
Automated Y-Bus Formation: Includes line charging (B/2) and transformer tap settings.
Q-Limit Checking: Includes logic to monitor reactive power limits at PV buses to ensure realistic generator operation.
Sparse Jacobian Matrix: Efficiently constructs the Jacobian matrix using partial derivatives of $P$ and $Q$.
Polar/Rectangular Utilities: Helper functions for coordinate conversions.

**🛠 Mathematical Overview**

The Newton-Raphson method solves the non-linear power flow equations by iteratively updating the voltage magnitudes and angles. The core of the algorithm is the Jacobian Matrix ($J$):$$\begin{bmatrix} \Delta P \\ \Delta Q \end{bmatrix} = \begin{bmatrix} J_1 & J_2 \\ J_3 & J_4 \end{bmatrix} \begin{bmatrix} \Delta \delta \\ \Delta |V| \end{bmatrix}$$Power Flow EquationsFor each bus $i$, the real ($P_i$) and reactive ($Q_i$) power injections are calculated as:$P_i = \sum_{k=1}^{n} |V_i||V_k|(G_{ik}\cos\theta_{ik} + B_{ik}\sin\theta_{ik})$$Q_i = \sum_{k=1}^{n} |V_i||V_k|(G_{ik}\sin\theta_{ik} - B_{ik}\cos\theta_{ik})$

📝 Dependencies
MATLAB (Tested on R2021a or later)

No additional toolboxes required.

🤝 Contributing
Feel free to fork this project and submit pull requests for enhancements, such as:

Adding Fast Decoupled Load Flow (FDLF) scripts.

Adding visualization for the power system graph.

Implementing line flow and loss calculations.
