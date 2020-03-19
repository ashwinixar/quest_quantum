/*
*****************************************************************************
*	About: The code implements Grover's Algorithm to solve 3-SAT problem	*
*	The	3-SAT function is same as in below link								*
*	https://qiskit.org/textbook/ch-applications/satisfiability-grover.html	*
*	Usage: Run in command prompt "3sat_grover.exe"							*
*	Author: Ashwini Kumar Malviya											*
*	Email: ashwinixar@gmail.com												*
*****************************************************************************
*/

//The 3-SAT function is f(x1, x2, x3) = (x1 ∨ x2 ∨ ¬x3) ∧ (¬x1 ∨ ¬x2 ∨ ¬x3) ∧ (¬x1 ∨ x2 ∨ x3)
//There exist t = 5 solutions for the above function f to be 1 i.e. f(x1, x2, x3) = 1
//Solutions are (x1,x2,x3) : (0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)
//Thus, Grover's iteration will be repeated 1 time i.e. times = (3.14 / 4) * sqrt(N / t) = (3.14 / 4) * sqrt(8 / 5) = 0.993

#include <stdio.h>
#include "QuEST.h"

int main (int narg, char *varg[])
{
    const int n = 3; //Input size (in bits) of the function "f"
    const int terms = 3; //Number of terms ANDed in the function "f"

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg((n + terms), env); //"n" qubits for function input and "terms" ancilla qubit
    initZeroState(qubits);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Create superpostion of x1, x2, and x3 on the first 3-qubits
    for(int i = 0; i < n; i++)
        hadamard(qubits, i);
    
    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

    //Function evaluation (term wise) and store the result in "terms" qubits
    int ctrls[n];
    for(int i = 0; i < n; i++)
    	ctrls[i] = i;
    //¬x implies 1 and x implies 0
    int term_form[][3] = { { 0, 0, 1 }, { 1, 1, 1 }, { 1, 0, 0 } };
    for(int i = 0; i < terms; i++)
    	multiStateControlledUnitary(qubits, ctrls, term_form[i], n, (i + n), ux);

    int targs[n + terms];
    for(int i = 0; i < (n + terms); i++)
    	targs[i] = i;
    int targs_terms[terms];
    for(int i = 0; i < terms; i++)
    	targs_terms[i] = (i + n);

    int times = 2; 
    for(int gi = 0; gi < times; gi++)
    {
    	//Marking Circuit
    	multiControlledPhaseFlip(qubits, targs_terms, terms);
    	
        //Diffusion Circuit
        for(int i = 0; i < n; i++)
            hadamard(qubits, i);
        for(int i = 0; i < n; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, (n + terms));
        for(int i = 0; i < n; i++)
            pauliX(qubits, i);
        for(int i = 0; i < n; i++)
            hadamard(qubits, i);
    }
	
    printf("\nCircuit output:\n");

    qreal prob;
    for(int i = 0; i < 8; i++)
    {
    	prob = getProbAmp(qubits, i);
    	printf("Probability amplitude of (x1, x2, x3) = (%d, %d, %d) : %f\n", (i & 1), ((i & (1 << 1)) >> 1), ((i & (1 << 2)) >> 2), prob);
    }

    destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return 0;
}