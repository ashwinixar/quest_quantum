/*
*************************************************************
*	About: The code implements Quantum algorithm to find 	*
*	hamming weight of a given number						*
*	Usage: Run in command prompt "hamming.exe <number>"		*
*	Example: hamming.exe 12
*	NOTE: Code can be adjusted to find hamming weight of a	*
*	superposition of several bit strings					* 
*	Author: Ashwini Kumar Malviya							*
*	Email: ashwinixar@gmail.com								*
*************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("\nUsage: hamming.exe <number>");
		printf("\nExample: hamming.exe 12\n");
		return 0;
	}

	int num = atoi(varg[1]);

    QuESTEnv env = createQuESTEnv();

    const int pi = 3.14;

    //"n" represents the block size
    const int n = 4;
    const int w_n = 4;

    Qureg qubits = createQureg((n + w_n), env);
    initZeroState(qubits);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //First n-qubits are set to represent "num"
    int temp = num;
    int idx = 0;
    while(temp)
    {
    	if(temp & 1) pauliX(qubits, idx);
    	idx++;
    	temp >>= 1;
    }
    //pauliX(qubits, 0);

    //Apply QFT on the last (w_n)-qubits to create superposition
    //for(int i = n; i < (w_n + n); i++)
    //    hadamard(qubits, i);
    for(int i = (n + w_n - 1); i >= n; i--)
    {
    	hadamard(qubits, i);
    	//printf("\nOn qubit %d hadamard is applied", i);
    	int k = 2;
    	for(int j = i - 1; j >= n; j--)
    	{
    		controlledPhaseShift(qubits, j, i, 2 * pi / (pow(2, k)));
    		//printf("\nOn qubit %d R%d is applied with control qubit %d", i, k, j);
    		k++;
    	}
    }

    for(int i = n - 1; i >= 0; i--)
    {
    	int k = 1;
    	for(int j = n; j < (w_n + n); j++)
    	{
    		controlledPhaseShift(qubits, i, j, 2 * pi / (pow(2, k)));
    		//printf("\nOn qubit %d R%d is applied with control qubit %d", j, k, i);
    		k++;
    	}
    }
    
    //Inverse QFT on the last (w_n)-qubits
    for(int i = n; i < (w_n + n); i++)
    {
    	for(int j = n; j < i; j++)
    	{
    		qreal angle = -2 * pi / (pow(2, i - j + 1));
    		controlledPhaseShift(qubits, j, i, angle);
    		//printf("\nOn qubit %d R%d(Inverse) is applied with control qubit %d", i, i - j + 1, j);
    	}
    	hadamard(qubits, i);
    	//printf("\nOn qubit %d hadamard is applied", i);
    }

    printf("\nCircuit output:\n");
    //Measures the qubits to get the value of last (w_n)-qubits
    qreal prob;
    int hamming_weight = 0;
    for(int i = n; i < (w_n + n); i++)
    {
        int outcome = measureWithStats(qubits, i, &prob);
        if(outcome) hamming_weight = hamming_weight ^ (0x1 << (i - n));
    	printf("Qubit %d collapsed to %d with probability %f\n", i, outcome, prob);
    }

    //Print the measured result
    printf("\nThe hamming weight of %d is %d\n", num, hamming_weight);

    destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return 0;
}
