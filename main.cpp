//Nicolas Stoian
//This program needs 3 command line arguments
//argv[1] "input1" for text file representing the list of points in x-y coordinates
//argv[2] "input2" for integer representing the number of K clusters
//argv[3] "output1" to write the 2-D image

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <limits>

using namespace std;

void loadPointSet(int** pointSet, ifstream& inFile);
void displayPointSet(int** pointSet, int numRows, int numCols);
void outputPointSet(int** pointSet, int numRows, int numCols, ofstream& outFile);
void outputPointSetWide(int** pointSet, int numRows, int numCols, ofstream& outFile);
void partitionPointSet(int** pointSet, int numRows, int numCols, int* numPoints, int k);
void computeCentroids(int** pointSet, int numRows, int numCols, int* numPoints, double* centroidX, double* centroidY, int k );
int checkDistances(int** pointSet, int numRows, int numCols, int* numPoints, double* centroidX, double* centroidY, int k );
double computeDistance(int x1, int y1, double x2, double y2);
int findMinI(double* distances, int k);

int main(int argc, char* argv[]){
    ifstream inFile;
    inFile.open(argv[1]);
    int numRows;
    inFile >> numRows;
    int numCols;
    inFile >> numCols;
    int** pointSet;
    pointSet = new int* [numRows];
    for(int i = 0; i < numRows; i++){
        pointSet[i] = new int [numCols];
    }
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            pointSet[row][col] = -1;
        }
    }
    loadPointSet(pointSet, inFile);
    stringstream numSets(argv[2]);
    int k;
    numSets >> k;
    int* numPoints = new int [k];
    double* centroidX = new double [k];
    double* centroidY = new double [k];
    for(int i = 0; i < k; i++){
        numPoints[i] = 0;
        centroidX[i] = 0;
        centroidY[i] = 0;
    }
    int numChangedLabel = 1;
    int iteration = 0;
    partitionPointSet(pointSet, numRows, numCols, numPoints, k);
    ofstream outFile;
    outFile.open(argv[3]);
    outFile << "iteration --> " << iteration << endl;
    outFile << "original partition of points, k = " << k << endl;
    outputPointSet(pointSet, numRows, numCols, outFile);
    while(numChangedLabel != 0){
        computeCentroids(pointSet, numRows, numCols, numPoints, centroidX, centroidY, k);
        numChangedLabel = checkDistances(pointSet, numRows, numCols, numPoints, centroidX, centroidY, k);
        iteration++;
        outFile << "iteration --> " << iteration << endl;
        outFile << "numChangedLabel --> " << numChangedLabel << endl;
        outputPointSet(pointSet, numRows, numCols, outFile);
    }
    inFile.close();
    outFile.close();
    for(int i = 0; i < numRows; i++){
        delete pointSet[i];
    }
    delete pointSet;
    delete numPoints;
    delete centroidX;
    delete centroidY;
}

void loadPointSet(int** pointSet, ifstream& inFile){
    int row;
    int col;
    while(inFile >> row){
        inFile >> col;
        pointSet[row][col] = 0;
    }
}

void displayPointSet(int** pointSet, int numRows, int numCols){
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] == -1){
                cout << " ";
            }
            else{
                cout << pointSet[row][col];
            }
        }
        cout << endl;
    }
}

void outputPointSet(int** pointSet, int numRows, int numCols, ofstream& outFile){
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] == -1){
                outFile << " ";
            }
            else{
                outFile << pointSet[row][col];
            }
        }
        outFile << endl;
    }
}

void partitionPointSet(int** pointSet, int numRows, int numCols, int* numPoints, int k){
    srand(time(NULL));
    int numPointsTotal = 0;
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] == 0){
                numPointsTotal++;
            }
        }
    }
    int* pointsLimit = new int [k];
    for(int i = 0; i < k; i++){
        pointsLimit[i] = numPointsTotal / k;
    }
    for(int i = 0; i < numPointsTotal % k; i++){
        pointsLimit[i]++;
    }
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] == 0){
                int group = rand() % k;
                while(numPoints[group] == pointsLimit[group]){
                    group = rand() % k;
                }
                pointSet[row][col] = group + 1;
                numPoints[group]++;
            }
        }
    }
    delete pointsLimit;
}

void computeCentroids(int** pointSet, int numRows, int numCols, int* numPoints, double* centroidX, double* centroidY, int k ){
    displayPointSet(pointSet, numRows, numCols);
    int* sumX = new int [k];
    int* sumY = new int [k];
    for(int i = 0; i < k; i++){
        sumX[i] = 0;
        sumY[i] = 0;
    }
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] > 0){
               sumX[pointSet[row][col] - 1] += row;
               sumY[pointSet[row][col] - 1] += col;
            }
        }
    }
    for(int i = 0; i < k; i++){
        centroidX[i] = 0;
        centroidY[i] = 0;
    }
    for(int i = 0; i < k; i++){
        if(numPoints[i] > 0){
            centroidX[i] = (double)sumX[i]/numPoints[i];
            centroidY[i] = (double)sumY[i]/numPoints[i];
        }
    }
    delete sumX;
    delete sumY;
}

int checkDistances(int** pointSet, int numRows, int numCols, int* numPoints, double* centroidX, double* centroidY, int k ){
    double* distances = new double [k];
    for(int i = 0; i < k; i++){
        distances[i] = 0;
    }
    int numChangedLabel = 0;
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(pointSet[row][col] > 0){
                for(int i = 0; i < k; i++){
                    if (centroidX[i] != 0 && centroidY[i] != 0){
                        distances[i] = computeDistance(row, col, centroidX[i], centroidY[i]);
                    }
                }
                int min_i = findMinI(distances, k);
                if(pointSet[row][col] != (min_i + 1)){
                    numPoints[pointSet[row][col] - 1]--;
                    numPoints[min_i]++;
                    pointSet[row][col] = min_i + 1;
                    numChangedLabel++;
                }
            }
        }
    }
    return numChangedLabel;
    delete distances;
}

double computeDistance(int x1, int y1, double x2, double y2){
    return sqrt(((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2)));
}

int findMinI(double* distances, int k){
    double minDistance = numeric_limits<double>::max();
    int min_i = 0;
    for(int i = 0; i < k; i++){
        if(distances[i] < minDistance && distances[i] != 0){
            minDistance = distances[i];
            min_i = i;
        }
    }
    return min_i;
}
