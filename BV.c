/*
*********************************************************************
*	About: The code implements BV Algorithm		 					*
*	Usage: Run in command prompt "BV.exe <hidden_string>"			*
*	<hidden_string> must be a 4-bit number (from 0 to 15)			*
*	Author: Ashwini Kumar Malviya									*
*	Email: ashwinixar@gmail.com										*
*********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include "QuEST.h"

//Dot product function
//Assumes that both x and y are of n-bits
int dot_product(int x, int y)
{
	int r = 0;
	while (x)
	{
		r ^= (x & 1) & (y & 1);
		x >>= 1;
		y >>= 1;
	}
	return r;
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("\nUsage: BV.exe <hidden_string>");
		printf("\n<hidden_string> must be a 4-bit number (from 0 to 15)\n");
		return 0;
	}

	//Let f : (F_2)^4 -> F_2 i.e. each input "x" to "f" is of 4-bits 
	//Let f(x) = s.x for some 4-bit string "s"
	int s = atoi(varg[1]); //Hidden string of 4-bits

	if(s > 15 || s < 0)
	{
		printf("\nUsage: BV.exe <hidden_string>");
		printf("\n<hidden_string> must be a 4-bit number (from 0 to 15)\n");
		return 0;
	}

	int oracle_f[16];
	for(int x = 0; x < 16; x++) oracle_f[x] = dot_product(x, s);

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(4, env);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    ComplexMatrixN e = createComplexMatrixN(4);
	for(int i = 0; i < 16; i++)
	{
	  	if(oracle_f[i] == 1) e.real[i][i] = -1;
	   	else e.real[i][i] = 1;
	}

	for(int i = 0; i < 4; i++)
        hadamard(qubits, i);

    int targs[] = { 0,1,2,3 };
    multiQubitUnitary(qubits, targs, 4, e);

    for(int i = 0; i < 4; i++)
        hadamard(qubits, i);

    qreal prob;
	int hs = 0;
	for(int i = 0; i < 4; i++)
	{
		int outcome = measureWithStats(qubits, i, &prob);
		if(outcome) hs ^= (outcome << i);
	}

	printf("\nThe hidden string is 0x%01X\n", hs);

	destroyQureg(qubits, env);
	destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

    return 0;
}