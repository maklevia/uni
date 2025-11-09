#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <iostream>
//#include "SerialFloyd.cpp"

int Min(int A, int B) {
    int Result = (A < B) ? A : B;
    if((A < 0) && (B >= 0)) Result = B;
    if((B < 0) && (A >= 0)) Result = A;
    if((A < 0) && (B < 0)) Result = -1;
    return Result;
}

// Function for the serial Floyd algorithm
void SerialFloyd(int *pMatrix, int Size) {
    int t1, t2;
    for(int k = 0; k < Size; k++)
        for(int i = 0; i < Size; i++)
            for(int j = 0; j < Size; j++)
                if((pMatrix[i * Size + k] != -1) &&
                (pMatrix[k * Size + j] != -1)) {
                    t1 = pMatrix[i * Size + j];
                    t2 = pMatrix[i * Size + k] + pMatrix[k * Size + j];
                    pMatrix[i * Size + j] = Min(t1, t2);

                }
}

// Function for comparing the matrices
bool CompareMatrices(int *pMatrix1, int *pMatrix2, int Size) {
    return std::equal(pMatrix1, pMatrix1 + Size * Size, pMatrix2);
}