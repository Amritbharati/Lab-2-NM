#include <iostream>
using namespace std;
int factorial(int n) {
    if (n <= 1)
        return 1;
    else
        return n * factorial(n - 1);
}

void calculateBackwardDifferenceTable(double x[], double y[], double backwardDifferenceTable[][10], int n) {
    for (int i = 0; i < n; i++)
        backwardDifferenceTable[i][0] = y[i];

    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            backwardDifferenceTable[i][j] = backwardDifferenceTable[i][j - 1] - backwardDifferenceTable[i - 1][j - 1];
        }
    }
}

double interpolate(double x[], double backwardDifferenceTable[][10], int n, double targetX) {
    double result = backwardDifferenceTable[n - 1][0];
    double term = (targetX - x[n - 1]) / (x[1] - x[0]);
    double productTerm = term;

    for (int i = 1; i < n; i++) {
        result += (productTerm * backwardDifferenceTable[n - 1][i]) / factorial(i);
        productTerm *= (term + i) / (i + 1);
    }

    return result;
}

int main() {
    int n;

    cout<< "Enter the number of data points: ";
    cin>> n;

    double x[10];
    double y[10];

    std::cout << "Enter the data points:\n";
    for (int i = 0; i < n; i++) {
        cout<< "x[" << i << "]: ";
        cin>> x[i];
        cout<< "y[" << i << "]: ";
        cin>> y[i];
    }

    double backwardDifferenceTable[10][10];

    calculateBackwardDifferenceTable(x, y, backwardDifferenceTable, n);

    double targetX;
    cout<< "Enter the value of x to interpolate: ";
    cin>> targetX;

    double interpolatedValue = interpolate(x, backwardDifferenceTable, n, targetX);

    cout<< "Interpolated value at x = " << targetX << ": " << interpolatedValue <<endl;

    return 0;
}


