Quantum Simulator in C
This repository contains a simple, command-line quantum computer simulator written in C. It provides a foundational understanding of quantum computing principles by simulating qubit state vectors and the effects of common quantum gates and algorithms.

🚀 Features
Qubit Simulation: Manages and updates the state of a multi-qubit system.

Core Gates: Includes implementations for fundamental quantum gates such as Hadamard (H), Pauli-X (X), Pauli-Y (Y), and Pauli-Z (Z).

Controlled-NOT (CNOT) Gate: Simulates the key two-qubit CNOT gate.

Probabilistic Measurement: Implements a function to measure the quantum state, collapsing the wave function to a single basis state.

💻 How to Compile and Run
To compile the simulator, you'll need a C compiler that supports C99 or later. The complex.h header is required for complex number operations. GCC is highly recommended.

Navigate to the project directory in your terminal.

Compile the source code using the following command:

gcc -o quantum_simulator main.c -lm

gcc: The C compiler.

-o quantum_simulator: Specifies the name of the output executable file.

main.c: The source file.

-lm: Links the math library, which is necessary for functions like sqrt() and M_PI.

Run the executable:

./quantum_simulator

⚛️ Implemented Quantum Algorithms
The simulator includes functions that demonstrate several classic quantum algorithms:

Deutsch–Jozsa Algorithm
This algorithm solves the Deutsch-Jozsa problem, which efficiently determines whether a given function is "constant" or "balanced." It accomplishes this with a single query to the quantum oracle, a task that would require many more queries classically.

Grover's Search Algorithm
Grover's algorithm is a quantum search algorithm that can find a specific item in an unstructured list with a quadratic speedup over classical algorithms. The simulator shows how the algorithm amplifies the amplitude of the target state, making it highly probable to measure.

Quantum Fourier Transform (QFT)
The QFT is a quantum analogue of the classical Discrete Fourier Transform. It is a crucial subroutine in many quantum algorithms, most notably Shor's algorithm for integer factorization. The simulator shows how the QFT transforms the state vector from the computational basis to the frequency basis.

Simon's Algorithm
Simon's algorithm solves the problem of finding a hidden period string, providing an exponential speedup over any classical algorithm. It was the first problem to demonstrate a clear computational advantage for quantum computers and laid the groundwork for Shor's algorithm.

🛠️ Code Structure
The main logic is contained within main.c, which includes key functions for:

initialize_state: Sets up a new quantum state vector.

apply_single_qubit_gate: Applies gates like H, X, Y, and Z.

apply_cnot_gate: Applies the CNOT gate.

apply_controlled_phase: A helper function for the QFT.

apply_qft: Implements the full QFT algorithm.

grovers_simulation: Runs a full simulation of Grover's algorithm.

deutsch_jozsa_simulation: Runs a full simulation of the Deutsch-Jozsa algorithm.

simons_simulation: Runs a full simulation of Simon's algorithm.

measure_state: Performs a single measurement.

print_state: Displays the current state vector.

➡️ Future Improvements
Entangled States: Implement functions to create and analyze entangled states like Bell pairs.

Additional Gates: Expand the gate library to include gates like Toffoli and Fredkin.

Modular Design: Separate the logic into multiple header and source files for better organization.

User Input: Add a feature for user-defined input for running gates and algorithms.
