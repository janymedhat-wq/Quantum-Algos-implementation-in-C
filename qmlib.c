#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

double complex multiply_complex(double complex a, double complex b)
{
    return a * b;
}

double complex add_complex(double complex a, double complex b)
{
    return a + b;
}

void initialize_state(int num_qubits, double complex *state_vector)
{
    int num_states = 1 << num_qubits;
    for (int i = 0; i < num_states; ++i)
    {
        state_vector[i] = 0.0 + 0.0 * I;
    }
    state_vector[0] = 1.0 + 0.0 * I;
}

void apply_single_qubit_gate(int num_qubits, double complex *state_vector,
                             double complex gate_matrix[2][2], int target_qubit)
{
    int num_states = 1 << num_qubits;
    int k = 1 << target_qubit;

    for (int i = 0; i < num_states; ++i)
    {
        if ((i & k) == 0)
        {
            int j = i | k;
            double complex temp0 = state_vector[i];
            double complex temp1 = state_vector[j];

            state_vector[i] = add_complex(multiply_complex(gate_matrix[0][0], temp0),
                                          multiply_complex(gate_matrix[0][1], temp1));
            state_vector[j] = add_complex(multiply_complex(gate_matrix[1][0], temp0),
                                          multiply_complex(gate_matrix[1][1], temp1));
        }
    }
}

void apply_cnot_gate(int num_qubits, double complex *state_vector,
                     int control_qubit, int target_qubit)
{
    int num_states = 1 << num_qubits;
    int control_mask = 1 << control_qubit;
    int target_mask = 1 << target_qubit;

    for (int i = 0; i < num_states; ++i)
    {
        if ((i & control_mask) == control_mask)
        {
            int j = i ^ target_mask;
            double complex temp = state_vector[i];
            state_vector[i] = state_vector[j];
            state_vector[j] = temp;
        }
    }
}

int measure_state(int num_qubits, const double complex *state_vector)
{
    int num_states = 1 << num_qubits;
    double r = (double)rand() / RAND_MAX;
    double cumulative_prob = 0.0;

    for (int i = 0; i < num_states; ++i)
    {
        double prob = cabs(state_vector[i]) * cabs(state_vector[i]);
        cumulative_prob += prob;
        if (r < cumulative_prob)
        {
            return i;
        }
    }
    return 0;
}

void print_state(int num_qubits, const double complex *state_vector)
{
    int num_states = 1 << num_qubits;
    for (int i = 0; i < num_states; ++i)
    {
        char binary_string[num_qubits + 1];
        for (int j = 0; j < num_qubits; ++j)
        {
            binary_string[num_qubits - 1 - j] = ((i >> j) & 1) ? '1' : '0';
        }
        binary_string[num_qubits] = '\0';
        printf("|%s> : (%.4f + %.4fi), Probability: %.4f\n",
               binary_string, creal(state_vector[i]), cimag(state_vector[i]),
               cabs(state_vector[i]) * cabs(state_vector[i]));
    }
    printf("\n");
}

void grovers_simulation(int num_qubits, int target_state)
{
    printf("--- Running Grover's Search for target state |%d> ---\n", target_state);
    int num_states = 1 << num_qubits;
    double complex *state_vector = (double complex *)malloc(num_states * sizeof(double complex));
    double complex H[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};

    initialize_state(num_qubits, state_vector);

    for (int i = 0; i < num_qubits; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
    }
    printf("Initial state after Hadamards:\n");
    print_state(num_qubits, state_vector);

    int num_iterations = (int)round(M_PI_4 * sqrt(num_states));
    printf("Number of iterations: %d\n\n", num_iterations);

    for (int iter = 0; iter < num_iterations; ++iter)
    {
        state_vector[target_state] *= -1;

        double complex mean = 0.0 + 0.0 * I;
        for (int i = 0; i < num_states; ++i)
        {
            mean = add_complex(mean, state_vector[i]);
        }
        mean /= (double)num_states;

        for (int i = 0; i < num_states; ++i)
        {
            state_vector[i] = add_complex(multiply_complex(2.0, mean), -state_vector[i]);
        }
    }

    printf("Final state after Grover's iterations:\n");
    print_state(num_qubits, state_vector);
    int result = measure_state(num_qubits, state_vector);
    printf("Grover's measurement result: %d\n", result);

    free(state_vector);
}

bool is_balanced(int input, int num_qubits)
{
    return (input >> (num_qubits - 2)) & 1;
}

void deutsch_jozsa_simulation(int num_qubits, bool func_is_constant)
{
    printf("--- Running Deutsch-Jozsa for a %s function ---\n", func_is_constant ? "constant" : "balanced");
    int num_states = 1 << num_qubits;
    double complex *state_vector = (double complex *)malloc(num_states * sizeof(double complex));
    double complex H[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};
    double complex X[2][2] = {{0, 1}, {1, 0}};

    initialize_state(num_qubits, state_vector);

    apply_single_qubit_gate(num_qubits, state_vector, X, num_qubits - 1);
    apply_single_qubit_gate(num_qubits, state_vector, H, num_qubits - 1);

    for (int i = 0; i < num_qubits - 1; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
    }
    printf("State after initial Hadamards:\n");
    print_state(num_qubits, state_vector);

    for (int i = 0; i < num_states; ++i)
    {
        bool output_bit;
        if (func_is_constant)
        {
            output_bit = true;
        }
        else
        {
            output_bit = (i >> (num_qubits - 2)) & 1;
        }

        if (output_bit)
        {
            int target_state_index = i;
            if ((target_state_index >> (num_qubits - 1)) & 1)
            {
                state_vector[target_state_index] *= -1;
            }
        }
    }
    printf("State after oracle application:\n");
    print_state(num_qubits, state_vector);

    for (int i = 0; i < num_qubits - 1; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
    }
    printf("Final state before measurement:\n");
    print_state(num_qubits, state_vector);

    int result = measure_state(num_qubits, state_vector);
    printf("Deutsch-Jozsa measurement result: %d\n", result);
    if (func_is_constant)
    {
        if (result == 0)
        {
            printf("Expected constant, got |0...0> as predicted.\n");
        }
        else
        {
            printf("Error in prediction.\n");
        }
    }
    else
    {
        if (result != 0)
        {
            printf("Expected balanced, got non-|0...0> as predicted.\n");
        }
        else
        {
            printf("Error in prediction.\n");
        }
    }

    free(state_vector);
}

void apply_controlled_phase(int num_qubits, double complex *state_vector,
                            int control_qubit, int target_qubit, double angle)
{
    int num_states = 1 << num_qubits;
    int control_mask = 1 << control_qubit;
    int target_mask = 1 << target_qubit;

    for (int i = 0; i < num_states; ++i)
    {
        if ((i & control_mask) == control_mask && (i & target_mask) == target_mask)
        {
            state_vector[i] = state_vector[i] * cexp(I * angle);
        }
    }
}

void apply_qft(int num_qubits, double complex *state_vector)
{
    double complex H[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};
    for (int i = 0; i < num_qubits; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
        for (int j = i + 1; j < num_qubits; ++j)
        {
            double angle = M_PI / (double)(1 << (j - i));
            apply_controlled_phase(num_qubits, state_vector, j, i, angle);
        }
    }

    for (int i = 0; i < num_qubits / 2; ++i)
    {
        double complex temp = state_vector[i];
        state_vector[i] = state_vector[num_qubits - 1 - i];
        state_vector[num_qubits - 1 - i] = temp;
    }
}

void simons_oracle(int num_qubits, double complex *state_vector, int s)
{
    int num_states = 1 << num_qubits;
    int input_qubits = num_qubits / 2;
    int output_qubits = num_qubits - input_qubits;
    double complex CNOT[2][2] = {{1, 0}, {0, 1}};
    double complex X[2][2] = {{0, 1}, {1, 0}};

    for (int i = 0; i < num_states; ++i)
    {
        int x = i >> output_qubits;
        int y_xor_x = x ^ s;
        int y = (i & ((1 << output_qubits) - 1)) ^ y_xor_x;
        int new_index = (x << output_qubits) | y;

        if (i < new_index)
        {
            double complex temp = state_vector[i];
            state_vector[i] = state_vector[new_index];
            state_vector[new_index] = temp;
        }
    }
}

void simons_simulation(int num_qubits, int s)
{
    printf("--- Running Simon's Algorithm for s = %d ---\n", s);
    int num_states = 1 << num_qubits;
    double complex *state_vector = (double complex *)malloc(num_states * sizeof(double complex));
    double complex H[2][2] = {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};

    initialize_state(num_qubits, state_vector);

    for (int i = 0; i < num_qubits; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
    }

    simons_oracle(num_qubits, state_vector, s);

    for (int i = 0; i < num_qubits / 2; ++i)
    {
        apply_single_qubit_gate(num_qubits, state_vector, H, i);
    }

    printf("Final state after Simon's algorithm:\n");
    print_state(num_qubits, state_vector);

    int result = measure_state(num_qubits, state_vector);
    printf("Simon's measurement result: %d\n", result);

    if ((result & s) == 0)
    {
        printf("Result is orthogonal to s, as expected.\n");
    }
    else
    {
        printf("Result is not orthogonal to s.\n");
    }

    free(state_vector);
}

int main()
{
    srand((unsigned int)time(NULL));

    double one_over_sqrt2 = 1.0 / sqrt(2.0);
    double complex H[2][2] = {{one_over_sqrt2, one_over_sqrt2}, {one_over_sqrt2, -one_over_sqrt2}};
    double complex X[2][2] = {{0, 1}, {1, 0}};
    double complex Y[2][2] = {{0, -1.0 * I}, {1.0 * I, 0}};
    double complex Z[2][2] = {{1, 0}, {0, -1}};

    printf("--- Single Qubit Simulation ---\n");
    int num_q_single = 1;
    double complex *q_state = (double complex *)malloc((1 << num_q_single) * sizeof(double complex));
    initialize_state(num_q_single, q_state);
    printf("Initial state of 1 qubit:\n");
    print_state(num_q_single, q_state);
    apply_single_qubit_gate(num_q_single, q_state, H, 0);
    printf("State after applying Hadamard gate:\n");
    print_state(num_q_single, q_state);
    int measured = measure_state(num_q_single, q_state);
    printf("Measured state: |%d>\n\n", measured);
    free(q_state);

    printf("--- Two Qubit CNOT Gate ---\n");
    int num_q_cnot = 2;
    double complex *cnot_state = (double complex *)malloc((1 << num_q_cnot) * sizeof(double complex));
    initialize_state(num_q_cnot, cnot_state);
    apply_single_qubit_gate(num_q_cnot, cnot_state, H, 0);
    printf("State of |00> after H on qubit 0:\n");
    print_state(num_q_cnot, cnot_state);
    apply_cnot_gate(num_q_cnot, cnot_state, 0, 1);
    printf("State after CNOT with control 0 and target 1:\n");
    print_state(num_q_cnot, cnot_state);
    free(cnot_state);

    printf("--- QFT Simulation on 3 Qubits ---\n");
    int num_qft_qubits = 3;
    double complex *qft_state = (double complex *)malloc((1 << num_qft_qubits) * sizeof(double complex));
    initialize_state(num_qft_qubits, qft_state);
    qft_state[2] = 1.0;
    printf("Initial state for QFT: |010>\n");
    print_state(num_qft_qubits, qft_state);
    apply_qft(num_qft_qubits, qft_state);
    printf("Final state after QFT:\n");
    print_state(num_qft_qubits, qft_state);
    free(qft_state);

    int num_qubits_algo = 4;

    printf("\n--- Quantum Algorithm Simulations ---\n");

    deutsch_jozsa_simulation(num_qubits_algo, true);
    deutsch_jozsa_simulation(num_qubits_algo, false);

    grovers_simulation(3, 5);

    int num_simons_qubits = 6;
    int s = 5;
    simons_simulation(num_simons_qubits, s);

    return 0;
}
