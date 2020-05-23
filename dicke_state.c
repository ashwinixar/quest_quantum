//0-th qubit represent the LSB

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

typedef unsigned int _uint;

void print_bin(_uint a, _uint n)
{
	for(int i = n - 1; i >= 0; i--)
		if(a & (1 << i)) printf("1");
		else printf("0");
}

_uint hw(_uint n)
{
	_uint w = 0;
	while(n)
	{
		if(n & 1) w++;
		n >>= 1;
	}
	return w;
}

int main (int narg, char *varg[])
{
	_uint n = 5;
	_uint k = 3;

	if(n < k)
	{
		printf("Value of \"k\" must be smaller than or equal to \"n\"!\n");
		return 0;
	}

	QuESTEnv env = createQuESTEnv();
	Qureg qubits = createQureg(n, env);
	initZeroState(qubits);

	for(int qb = 0; qb < k; qb++)
		pauliX(qubits, qb);

	_uint _n = n;
	_uint base_qb = 0;
	for(_uint _k = k; _k > 0; )
	{
		_uint a_qb = base_qb + 1;
		
		int *ctrls = (int *)malloc(sizeof(int) * 2);
		int ctrls_size = 0;
		ctrls[ctrls_size++] = base_qb;
		ctrls[ctrls_size] = base_qb;
		int num = 1;
		for(_uint i = base_qb + 1; i <= base_qb + _k; i++)
		{
			double theta = 2.0 * acos(sqrt((double)num / _n));
			Complex alpha, beta;
			Vector axis = {0, 1, 0}; //y-axis for Rotation about y-axis
			alpha.real = cos(theta / 2.0);
			alpha.imag = -sin(theta / 2.0) * axis.z;
			beta.real = sin(theta / 2.0) * axis.y;
			beta.imag = -sin(theta / 2.0) * axis.x;

			ComplexMatrix2 U_Ry = {
				.real={{alpha.real,-beta.real},{beta.real,alpha.real}},
				.imag={{alpha.imag,beta.imag},{beta.imag,-alpha.imag}}
			};
			
			controlledNot(qubits, i, base_qb);
			multiControlledUnitary(qubits, ctrls, ctrls_size, a_qb, U_Ry);
			controlledNot(qubits, i, base_qb);
			
			a_qb++;
			ctrls_size = 2;
			ctrls[1]++;
			num++;
		}
		_n--;
		if(_n == _k) _k--;
		base_qb++;
	}

	qreal prob = 0.0;
	for(int i = 0; i < (int)pow(2, n); i++)
	{
		prob = getProbAmp(qubits, i);
		//if(prob == 0.0) continue;
		printf("Prob. of ");
		print_bin(i, n);
		printf(" (wt=%d) state is %f\n", hw(i), prob);
	}
	
	destroyQureg(qubits, env);
	destroyQuESTEnv(env);

	return 0;
}