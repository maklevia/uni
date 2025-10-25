#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <iostream>
#include <math.h>

int ProcNum = 0; // Number of the available processes
int ProcRank = 0; // Rank of the current process
int *pParallelPivotPos; // Number of rows selected as the pivot ones
int *pProcPivotIter; // Number of iterations, at which the process
int* pProcInd; // Number of the first row located on the processes
int* pProcNum; // Number of the linear system rows located on the processes

// rows were used as the pivot ones

// Function for simple initialization of the matrix and the vector elements
void DummyDataInitialization (double* pMatrix, double* pVector, int Size) {
    int i, j; // Loop variables
    for (i=0; i<Size; i++) {
        pVector[i] = i+1;
        for (j=0; j<Size; j++) {
            if (j <= i)
                pMatrix[i*Size+j] = 1;
            else
                pMatrix[i*Size+j] = 0;
        }
    }
}

// Function for random initialization of the matrix and the vector elements
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
    int i, j; // Loop variables
    srand(unsigned(clock()));
    for (i=0; i<Size; i++) {
        pVector[i] = rand()/double(1000);
        for (j=0; j<Size; j++) {
            if (j <= i)
                pMatrix[i*Size+j] = rand()/double(1000);
            else
                pMatrix[i*Size+j] = 0;
        }
    }
}
void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
{
    int i, j; // Loop variables
    for (i = 0; i < RowCount; i++)
    {
        for (j = 0; j < ColCount; j++)
            printf("%7.4f ", pMatrix[i * ColCount + j]);
        printf("\n");
    }
}
// Function for formatted vector output
void PrintVector(double *pVector, int Size)
{
    int i;
    for (i = 0; i < Size; i++)
        printf("%7.4f ", pVector[i]);
}

// Function for memory allocation and data initialization
void ProcessInitialization (double* &pMatrix, double* &pVector, double* &pResult, double* &pProcRows, double* &pProcVector,
double* &pProcResult, int &Size, int &RowNum) {

    int RestRows = Size;
    for (int i=0; i<ProcRank; i++)
        RestRows = RestRows-RestRows/(ProcNum-i);
    RowNum = RestRows/(ProcNum-ProcRank);
    pProcRows = new double [RowNum*Size];
    pProcVector = new double [RowNum];
    pProcResult = new double [RowNum];
    if (ProcRank == 0) {
        pMatrix = new double [Size*Size];
        pVector = new double [Size];
        pResult = new double [Size];
        // Initialization of the matrix and the vector elements
        //DummyDataInitialization (pMatrix, pVector, Size);
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}

// Function for computational process termination
void ProcessTermination (double* pMatrix, double* pVector, double* pResult,
double* pProcRows, double* pProcVector, double* pProcResult) {
    if (ProcRank == 0) {
        delete [] pMatrix;
        delete [] pVector;
        delete [] pResult;
    }
    delete [] pProcRows;
    delete [] pProcVector;
    delete [] pProcResult;
}

// Function for the data distribution among the processes
void DataDistribution(double* pMatrix, double* pProcRows, double* pVector,
                      double* pProcVector, int Size, int RowNum) {
    int *pSendNumElems; // Number of elements (for matrix) sent to the process
    int *pSendIndElems; // Index (in elements) of the first matrix element sent
    int *pSendNumRows;  // Number of rows (for vector) sent to the process
    int *pSendIndRows;  // Index (in rows) of the first row sent
    int RestRows = Size;
    int i;

    // Allocate arrays
    pSendIndElems = new int [ProcNum];
    pSendNumElems = new int [ProcNum];
    pSendIndRows = new int [ProcNum];
    pSendNumRows = new int [ProcNum];
    pProcInd = new int [ProcNum]; // first row index (in rows) on each process
    pProcNum = new int [ProcNum]; // number of rows on each process

    // Compute row distribution (balanced progressive algorithm)
    RestRows = Size;
    int startRow = 0;
    for (i = 0; i < ProcNum; i++) {
        int rowsForProc = RestRows / (ProcNum - i);
        pSendNumRows[i] = rowsForProc;         // rows (for vector)
        pSendIndRows[i] = startRow;           // start row index
        pProcInd[i] = pSendIndRows[i];        // store as rows
        pProcNum[i] = pSendNumRows[i];        // store as rows
        // for matrix we need counts in elements
        pSendNumElems[i] = pSendNumRows[i] * Size;
        pSendIndElems[i] = pSendIndRows[i] * Size;
        startRow += rowsForProc;
        RestRows -= rowsForProc;
    }

    // Scatter the matrix rows (counts are in elements)
    MPI_Scatterv(pMatrix, pSendNumElems, pSendIndElems, MPI_DOUBLE,
                 pProcRows, pSendNumElems[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Scatter the vector blocks (counts are in rows)
    MPI_Scatterv(pVector, pSendNumRows, pSendIndRows, MPI_DOUBLE,
                 pProcVector, pSendNumRows[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Free temporary arrays for matrix/vector displacements
    delete [] pSendNumElems;
    delete [] pSendIndElems;
    delete [] pSendNumRows;
    delete [] pSendIndRows;
}



// Function for the column elimination
void ParallelEliminateColumns(double* pProcRows, double* pProcVector,
double* pPivotRow, int Size, int RowNum, int Iter) {
    double PivotFactor;
    for (int i=0; i<RowNum; i++) {
        if (pProcPivotIter[i] == -1) {
            PivotFactor = pProcRows[i*Size+Iter] / pPivotRow[Iter];
            for (int j=Iter; j<Size; j++) {
                pProcRows[i*Size + j] -= PivotFactor* pPivotRow[j];
            }
            pProcVector[i] -= PivotFactor * pPivotRow[Size];
        }
    }
}

// Function for the Gaussian elimination
// Function for the Gaussian elimination
void ParallelGaussianElimination (double* pProcRows, double* pProcVector,
int Size, int RowNum) {
    double MaxValue; // Value of the pivot element of thе process
    int PivotPos; // Position of the pivot row in the process stripe
    struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;
    double *pPivotRow; // Pivot row of the current iteration
    pPivotRow = new double [Size+1];

    // The iterations of the Gaussian elimination
    for (int i=0; i<Size; i++) {
        MaxValue = 0.0;
        // Calculating the local pivot row
        for (int j=0; j<RowNum; j++) {
            if ((pProcPivotIter[j] == -1) &&
            (MaxValue < fabs(pProcRows[j*Size+i]))) {
                MaxValue = fabs(pProcRows[j*Size+i]);
                PivotPos = j;
            }
        }
        // Finding the global pivot row
        ProcPivot.MaxValue = MaxValue;
        ProcPivot.ProcRank = ProcRank;
        // Finding the pivot process
        MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
        MPI_COMM_WORLD);

        // Storing the number of the pivot row
        if ( ProcRank == Pivot.ProcRank ){
            pProcPivotIter[PivotPos]= i;
            pParallelPivotPos[i]= pProcInd[ProcRank] + PivotPos;
        }
        MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank,
        MPI_COMM_WORLD);

        // Broadcasting the pivot row
        if ( ProcRank == Pivot.ProcRank ){
            // Fill the pivot row
            for (int j=0; j<Size; j++) {
                pPivotRow[j] = pProcRows[PivotPos*Size + j];
            }
            pPivotRow[Size] = pProcVector[PivotPos];
        }
        MPI_Bcast(pPivotRow, Size+1, MPI_DOUBLE, Pivot.ProcRank,
        MPI_COMM_WORLD);
        // Column elimination
        ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow, Size,
        RowNum, i);
    }
    delete [] pPivotRow;
}

// Function for finding the pivot row of the back substitution
void FindBackPivotRow(int RowIndex, int &IterProcRank,
int &IterPivotPos) {
    for (int i=0; i<ProcNum-1; i++) {
        if ((pProcInd[i]<=RowIndex) && (RowIndex<pProcInd[i+1]))
            IterProcRank = i;
    }
    if (RowIndex >= pProcInd[ProcNum-1])
        IterProcRank = ProcNum-1;
    IterPivotPos = RowIndex - pProcInd[IterProcRank];
}

// Function for the back substitution
void ParallelBackSubstitution (double* pProcRows, double* pProcVector,
double* pProcResult, int Size, int RowNum) {
    int IterProcRank; // Rank of the process with the current pivot row
    int IterPivotPos; // Position of the pivot row of the process
    double IterResult; // Calculated value of the current unknown
    double val;
    // The iterations of the back substitution
    for (int i=Size-1; i>=0; i--) {
        // Calculating the rank of the process, which holds the pivot row
        FindBackPivotRow(pParallelPivotPos[i],IterProcRank, IterPivotPos);
        // Calculating the unknown
        if (ProcRank == IterProcRank) {
            IterResult = pProcVector[IterPivotPos] /
            pProcRows[IterPivotPos*Size+i];
            pProcResult[IterPivotPos] = IterResult;
        }
        // Broadcasting the value of the current unknown
        MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);
        // Updating the values of the vector
        for (int j=0; j<RowNum; j++)
            if ( pProcPivotIter[j] < i ) {
                val = pProcRows[j*Size + i] * IterResult;
                pProcVector[j]=pProcVector[j] - val;
            }
    }
}

// Function for execution of the parallel Gauss algorithm
void ParallelResultCalculation(double* pProcRows, double* pProcVector,
double* pProcResult, int Size, int RowNum) {
    // Memory allocation
    pParallelPivotPos = new int [Size];
    pProcPivotIter = new int [RowNum];
    for (int i=0; i<RowNum; i++)
        pProcPivotIter[i] = -1;

    // Gaussian elimination
    ParallelGaussianElimination (pProcRows, pProcVector, Size, RowNum);
    // Back substitution
    ParallelBackSubstitution (pProcRows, pProcVector, pProcResult, Size, RowNum);

    // Memory deallocation
   // delete [] pParallelPivotPos;
    //delete [] pProcPivotIter;
}

// Function for gathering the result vector
void ResultCollection(double* pProcResult, double* pResult) {
    //Gathering the result vector on the pivot processor
    MPI_Gatherv(pProcResult, pProcNum[ProcRank], MPI_DOUBLE, pResult,
    pProcNum, pProcInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Function for testing the data distribution
void TestDistribution(double* pMatrix, double* pVector, double* pProcRows,
double* pProcVector, int Size, int RowNum) {
    if (ProcRank == 0) {
        printf("Initial Matrix: \n");
        PrintMatrix(pMatrix, Size, Size);
        printf("Initial Vector: \n");
        PrintVector(pVector, Size);
    }
    for (int i=0; i<ProcNum; i++) {
        if (ProcRank == i) {
            printf("\nProcRank = %d \n", ProcRank);
            printf(" Matrix Stripe:\n");
            PrintMatrix(pProcRows, RowNum, Size);
            printf(" Vector: \n");
            PrintVector(pProcVector, RowNum);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Function for testing the result
void TestResult(double* pMatrix, double* pVector, double* pResult,
int Size) {
    /* Buffer for storing the vector, that is a result of multiplication
    of the linear system matrix by the vector of unknowns */
    double* pRightPartVector;
    // Flag, that shows wheather the right parts vectors are identical or not
    int equal = 0;
    double Accuracy = 1.e-6; // Comparison accuracy
    if (ProcRank == 0) {
        pRightPartVector = new double [Size];
        for (int i=0; i<Size; i++) {
            pRightPartVector[i] = 0;
            for (int j=0; j<Size; j++) {
                pRightPartVector[i] +=
                pMatrix[i*Size+j]*pResult[pParallelPivotPos[j]];
            }
        }
        for (int i=0; i<Size; i++) {
            if (fabs(pRightPartVector[i]-pVector[i]) > Accuracy)
            equal = 1;
        }
        if (equal == 1)
            printf("The result of the parallel Gauss algorithm is NOT correct."
            "Check your code.");
        else
            printf("The result of the parallel Gauss algorithm is correct.");
        delete [] pRightPartVector;
    }
}

// Function for formatted result vector output
void PrintResultVector (double* pResult, int Size) {
    int i;
    for (i=0; i<Size; i++)
        printf("%7.4f ", pResult[pParallelPivotPos[i]]);
}


int main(int argc, char *argv[]) {
    double* pMatrix; // Matrix of the linear system
    double* pVector; // Right parts of the linear system
    double* pResult; // Result vector
    double *pProcRows; // Rows of the matrix A
    double *pProcVector; // Block of the vector b
    double *pProcResult; // Block of the vector x
    int Size; // Size of the matrix and the vectors
    int RowNum; // Number of the matrix rows
    double Start, Finish, Duration;

    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    const int TEST_SIZES[] = {10, 100, 500, 1000, 1500, 2000, 2500, 3000};
    const int NUM_TESTS = sizeof(TEST_SIZES) / sizeof(TEST_SIZES[0]);

    if (ProcRank == 0)
        printf("Parallel Gauss algorithm for solving linear systems\n");
    for (int test = 0; test < NUM_TESTS; ++test) {
        Size = TEST_SIZES[test];
        if (ProcRank == 0)
            printf("\n=== Running test with matrix size %d ===\n", Size);

        // Memory allocation and data initialization
        ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcVector,
        pProcResult, Size, RowNum);

        Start = MPI_Wtime();
        // Distributing the initial data between the processes
        DataDistribution(pMatrix, pProcRows, pVector, pProcVector, Size, RowNum);
        //TestDistribution(pMatrix, pVector, pProcRows, pProcVector, Size, RowNum);

        // The execution of the parallel Gauss algorithm
        ParallelResultCalculation (pProcRows, pProcVector, pProcResult, Size, RowNum);
        //TestDistribution(pMatrix, pVector, pProcRows, pProcVector, Size, RowNum);

        // Gathering the result vector
        ResultCollection(pProcResult, pResult);
        //if (ProcRank == 0) {
        // printf ("Result vector \n");
        //PrintResultVector(pResult, Size);
        // }
        Finish = MPI_Wtime();
        Duration = Finish-Start;

        // Testing the result
        //TestResult(pMatrix, pVector, pResult, Size);
        // Printing the time spent by parallel Gauss algorithm
        if (ProcRank == 0)
            printf("\n Time of execution: %f\n", Duration);

        // Process termination
        ProcessTermination (pMatrix, pVector, pResult, pProcRows, pProcVector,
        pProcResult);
    }
    MPI_Finalize();
    return 0;
}

