#include <stdio.h>
#include <stdlib.h>
// #include <conio.h>
#include <time.h>
// Function for simple setting the matrix and vector elements
void DummyDataInitialization (double* pMatrix, double* pVector, int Size) {
    int i, j; // Loop variables
    for (i=0; i<Size; i++) {
        pVector[i] = 1;
        for (j=0; j<Size; j++)
            pMatrix[i*Size+j] = i;
    }
}
// Function for random setting the matrix and vector elements
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
    int i, j; // Loop variables
    srand(unsigned(clock()));
    for (i=0; i<Size; i++) {
        pVector[i] = rand()/double(1000);
        for (j=0; j<Size; j++)
            pMatrix[i*Size+j] = rand()/double(1000);
    }
}
// Function for memory allocation and data initialization
void ProcessInitialization (double* &pMatrix, double* &pVector,
double* &pResult, int &Size) {
    // Setting the size of the initial matrix and vector
    printf("\nEnter the size of the initial objects: ");
    scanf("%d", &Size);
    if (Size <= 0) {
        return;
    }
    printf("\nChosen objects size = %d\n", Size);

    // Memory allocation
    pMatrix = new double [Size*Size];
    pVector = new double [Size];
    pResult = new double [Size];
    // Setting the values of the matrix and vector elements
    //DummyDataInitialization(pMatrix, pVector, Size);

    // Random definition of matrix and vector elements
    RandomDataInitialization(pMatrix, pVector, Size);
}
// Function for formatted matrix output
void PrintMatrix (double* pMatrix, int RowCount, int ColCount) {
    int i, j; // Loop variables
    for (i=0; i<RowCount; i++) {
        for (j=0; j<ColCount; j++)
            printf("%7.4f ", pMatrix[i*RowCount+j]);
        printf("\n");
    }
}
void PrintVector (double* pVector, int Size) {
    int i;
    for (i=0; i<Size; i++)
        printf("%7.4f ", pVector[i]);
}
// Function for matrix-vector multiplication
void ResultCalculation(double* pMatrix, double* pVector, double* pResult,
int Size) {
    int i, j; // Loop variables
    for (i=0; i<Size; i++) {
        pResult[i] = 0;
        for (j=0; j<Size; j++)
            pResult[i] += pMatrix[i*Size+j]*pVector[j];
    }
}
// Function for computational process termination
void ProcessTermination(double* pMatrix,double* pVector,double* pResult) {
    delete [] pMatrix;
    delete [] pVector;
    delete [] pResult;
}
int main() {
    double* pMatrix; // First argument - initial matrix
    double* pVector; // Second argument - initial vector
    double* pResult; // Result vector for matrix-vector multiplication
    int Size; // Sizes of initial matrix and vector
    time_t start, finish;
    double duration;
    printf("Serial matrix-vector multiplication program\n");

    // Memory allocation and data initialization
    ProcessInitialization(pMatrix, pVector, pResult, Size);
    if (Size <= 0) {
        printf("\nSize of objects must be greater than 0!\n");
        return 0;
    }


    // Matrix and vector output
    // printf ("Initial Matrix \n");
    // PrintMatrix(pMatrix, Size, Size);
    // printf("Initial Vector \n");
    // PrintVector(pVector, Size);
    // Matrix-vector multiplication

    start = clock();
    ResultCalculation(pMatrix, pVector, pResult, Size);
    finish = clock();
    duration = (finish-start)/double(CLOCKS_PER_SEC);
    // Printing the result vector
    // printf ("\n Result Vector: \n");
    // PrintVector(pResult, Size);

    // Printing the time spent by matrix-vector multiplication
    printf("\n Time of execution: %f\n", duration);
    // Computational process termination
    ProcessTermination(pMatrix, pVector, pResult);
}