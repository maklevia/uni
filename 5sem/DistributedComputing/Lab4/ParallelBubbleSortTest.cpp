#include <algorithm>
#include <cstdio>
using namespace std;

// Function for the serial bubble sort algorithm
void SerialBubbleSort(double *pData, int DataSize) {
    double Tmp;
    for(int i = 1; i < DataSize; i++)
        for(int j = 0; j < DataSize - i; j++)
            if(pData[j] > pData[j + 1]) {
                Tmp = pData[j];
                pData[j] = pData[j + 1];
                pData[j + 1] = Tmp;
            }
}

// Function for copying the sorted data
void CopyData(double *pData, int DataSize, double *pDataCopy) {
copy(pData, pData + DataSize, pDataCopy);
}

// Function for comparing the data
bool CompareData(double *pData1, double *pData2, int DataSize) {
return equal(pData1, pData1 + DataSize, pData2);
}