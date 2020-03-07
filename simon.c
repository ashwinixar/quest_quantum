/*
*****************************************************
*	About: The code implements Simon's Algorithm	*
*	Usage: Run in command prompt "simon.exe"		*
*	Author: Ashwini Kumar Malviya					*
*	Email: ashwinixar@gmail.com						*
*****************************************************
*/

//Let the function f : (F_2)^3 -> (F_2)^3, bits are indexed from left to right
//Each component function f_i of f is defined as:
//f_1 = x1 XOR x2 XOR x3 XOR 1
//f_2 = x3
//f_3 = x1x3 XOR x2x3 XOR x1 XOR x2 XOR x3 XOR 1
//f(000) = 101, f(001) = 010, f(010) = 000, f(011) = 110, f(100) = 000, f(101) = 110, f(110) = 101, f(111) = 010
//Below array declaration not needed, can be commented after commenting code at lines 164-165, also uncomment the code at lines 154-162
int f[] = { 5, 2, 0, 6, 0, 6, 5, 2 };

//NOTE: Since the function is small i.e. only 3-bits wide (3-bit input and output), if you are lucky then you may find the solution ..
// .. or the code may stuck in a loop or it will return no solution found. Break the code and rerun it to find the solution.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"


//Dot product function
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

//Before insertion checks if the vector is already present or not and check linear independency
int is_in(int *x, int n, int a)
{
	for(int i = 0; i < n; i++)
		if(x[i] == a && dot_product(x[i], a) == 0) return 1;
	return 0;
}

//Given "A" and "b", the function solves "Ax = b" for all non-trivial "x" (i.e. x != 0)
//"n" elements (or equations are present in "A[]")
//"width" represent that there are "width" number of variable
//"bit = { 0, 1 }", solve for "b" where each element in it is set to "bit"
//Find all the non-trivial vectors "x" and store in "solutions"
//Return "count" as number of solutions found
int binary_lin_system(int A[], int n, int width, int bit, int **solutions)
{
	int count = 0;
	*solutions = (int *)calloc((int)pow(2, width), sizeof(int));
	for (int x = 1; x < (int)pow(2, width); x++)
	{
		int found = 1;
		for (int i = 0; i < n; i++)
		{
			if (dot_product(A[i], x) != bit)
			{
				found = 0;
				break;
			}
		}
		if (found) (*solutions)[count++] = x;
	}
	return count;
}


int main (int narg, char *varg[])
{
    int n = 3; //Input and Output size (in bits) of the function "f"

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(n + n, env);
    initZeroState(qubits);
    reportQuregParams(qubits);
    reportQuESTEnv(env);
    
    //Stores (n - 1) linearly independent n-bit strings
    int *vec = (int *)calloc(n, sizeof(int));
    int vec_size = 0;

    //2-qubit Controlled-Not on the 3rd qubit
    ComplexMatrixN u = createComplexMatrixN(3);
	u.real[6][7] = 1;
	u.real[7][6] = 1;
	for (int i = 0; i < 6; i++)
	   	u.real[i][i] = 1;

    //Repeat until (n - 1) linearly independent n-bit strings are found
	//The function is only 3-bits wide, so in few iterations we get the answer "s", thus run the loop until first non-zero solution is found
    while(vec_size < 1) //For wider functions change the condition to (vec_size <= (n - 1))
    {
    	//STARTS - Simon's Algorithm

	    for(int i = 0; i < n; i++)
		    hadamard(qubits, i);

		//Function evaluation circuit
		//Qubits indexed from (n - 1) to 0 represents x1, x2, ..., xn
		//Qubits indexed from (2n - 1) to n represents the function output bits f1, f2, ..., fn
		//"x1", "x2", and "x3" stores the qubit index representing x1, x2, and x3, respectively i.e. 0th qubit represent x1, 1th index qubit represents x2, and so on.
		int x1 = 0, x2 = 1, x3 = 2;
		//"f_bit" stores the bit position for corresponding f_i
		int f_bit;
		//Evaluation of f1
		f_bit = (n + n - 1);
		pauliX(qubits, f_bit);
		controlledNot(qubits, x1, f_bit);
		controlledNot(qubits, x2, f_bit);
		controlledNot(qubits, x3, f_bit);
		//Evaluation of f2
		f_bit = (n + n - 2);
		controlledNot(qubits, x3, f_bit);
		//Evaluation of f3
		f_bit = n; //(n + n - 3) = m
		pauliX(qubits, f_bit);
		controlledNot(qubits, x1, f_bit);
		controlledNot(qubits, x2, f_bit);
		controlledNot(qubits, x3, f_bit);
		
		int ctrls[] = { x1, x3, f_bit }; //Applying x1x3
	    multiQubitUnitary(qubits, ctrls, 3, u);
	    ctrls[0] = x2; //Applying x2x3
	    ctrls[1] = x3;
	    multiQubitUnitary(qubits, ctrls, 3, u);

		for(int i = 0; i < n; i++)
		    hadamard(qubits, i);
		
		//ENDS - Simon's Algorithm

		//Measuring all the qubits
		int outcome, y = 0;
		qreal prob;
		for(int i = 0; i < n; i++)
		{
			outcome = measureWithStats(qubits, i, &prob);
			if(outcome) y ^= (outcome << (n - 1 - i));
		}

		if(!is_in(vec, vec_size, y) && y != 0) vec[vec_size++] = y;
    }

    //Uncomment for wider functions
    /*
    printf("\n");
    int *s, s_size;
    s_size = binary_lin_system(vec, vec_size, n, 0, &s);
    for(int i = 0; i < s_size; i++)
    {
    	printf("%X is the hidden string\n", s[i]);
    }
    */
    //Comment for wider functions
    if(f[0] == f[vec[0]]) printf("\nThe hidden string is %X (hex)\n", vec[0]);
	else printf("\nNo solution found!\n");


	destroyComplexMatrixN(u);
	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return 0;
}