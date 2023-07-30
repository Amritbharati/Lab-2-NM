#include <iostream>
using namespace std;
int factorial(int n) {
    if (n <= 1)
        return 1;
    else
        return n * factorial(n - 1);
}
void calculateCentralDifferenceTable(double x[], double y[], double centralDifferenceTable[][10], int n) {
    for (int i = 0; i < n; i++)
        centralDifferenceTable[i][0] = y[i];

    for (int j = 1; j < n; j++) {
        for (int i = j; i < n; i++) {
            centralDifferenceTable[i][j] = centralDifferenceTable[i][j - 1] - centralDifferenceTable[i - 1][j - 1];
        }
    }
}
double interpolate(double x[], double centralDifferenceTable[][10], int n, double targetX) {
    double h = x[1] - x[0];
    int middleIndex = n / 2;
    double term = (targetX - x[middleIndex]) / h;
    double result = centralDifferenceTable[middleIndex][0];

    double productTerm = term;
    double denominator = 1;

    for (int i = 1; i < n; i++) {
        if (i % 2 == 0)
            denominator *= i;

        if (i % 2 != 0)
            result += (productTerm * centralDifferenceTable[middleIndex][(i - 1) / 2]) / denominator;
        else
            result += (productTerm * centralDifferenceTable[middleIndex][i / 2]) / denominator;

        productTerm *= (term * term - i * i) / ((i + 1) * h * h);
    }

    return result;
}

int main() {
    int n;

    cout<< "Enter the number of data points: ";
    cin>> n;

    double x[10];
    double y[10];

    cout<< "Enter the data points:\n";
    for (int i = 0; i < n; i++) {
        cout<< "x[" << i << "]: ";
        cin>> x[i];
        cout<< "y[" << i << "]: ";
        cin>> y[i];
    }

    double centralDifferenceTable[10][10];

    calculateCentralDifferenceTable(x, y, centralDifferenceTable, n);

    double targetX;
    cout<< "Enter the value of x to interpolate: ";
    cin>> targetX;

    double interpolatedValue = interpolate(x, centralDifferenceTable, n, targetX);

    cout<< "Interpolated value at x = " << targetX << ": " << interpolatedValue <<endl;

    return 0;
}


