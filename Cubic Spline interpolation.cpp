#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double x;
    double y;
} DataPoint;

typedef struct {
    double a, b, c, d;
} CubicCoefficients;

// Function to perform cubic spline interpolation
void cubicSplineInterpolation(DataPoint* data, int n, CubicCoefficients* coefficients) {
    double* h = (double*)malloc(sizeof(double) * n);
    double* alpha = (double*)malloc(sizeof(double) * n);
    double* l = (double*)malloc(sizeof(double) * (n + 1));
    double* mu = (double*)malloc(sizeof(double) * (n + 1));
    double* z = (double*)malloc(sizeof(double) * (n + 1));

    for (int i = 1; i < n; i++)
        h[i] = data[i].x - data[i - 1].x;

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n - 1; i++) {
        alpha[i] = 3.0 * (data[i + 1].y - data[i].y) / h[i] - 3.0 * (data[i].y - data[i - 1].y) / h[i - 1];
        l[i] = 2.0 * (data[i + 1].x - data[i - 1].x) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = 0;
    coefficients[n - 1].c = 0;

    for (int j = n - 2; j >= 0; j--) {
        coefficients[j].c = z[j] - mu[j] * coefficients[j + 1].c;
        coefficients[j].b = (data[j + 1].y - data[j].y) / h[j] - h[j] * (coefficients[j + 1].c + 2.0 * coefficients[j].c) / 3.0;
        coefficients[j].d = (coefficients[j + 1].c - coefficients[j].c) / (3.0 * h[j]);
        coefficients[j].a = data[j].y;
    }

    free(h);
    free(alpha);
    free(l);
    free(mu);
    free(z);
}

// Function to evaluate the cubic spline at a given x
double evaluateCubicSpline(double x, DataPoint* data, CubicCoefficients* coefficients, int n) {
    int i = 0;
    while (i < n && x > data[i].x)
        i++;

    i--; // The index of the interval containing x

    double dx = x - data[i].x;
    double result = coefficients[i].a + coefficients[i].b * dx + coefficients[i].c * dx * dx + coefficients[i].d * dx * dx * dx;
    return result;
}

int main() {
    int n;
    printf("Enter the number of data points: ");
    scanf("%d", &n);

    DataPoint* data = (DataPoint*)malloc(sizeof(DataPoint) * n);
    CubicCoefficients* coefficients = (CubicCoefficients*)malloc(sizeof(CubicCoefficients) * n);

    printf("Enter the data points (x, y):\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf %lf", &data[i].x, &data[i].y);
    }

    // Calculate cubic spline coefficients
    cubicSplineInterpolation(data, n, coefficients);

    // Evaluate the cubic spline at a given x
    double x;
    printf("Enter the value of x for interpolation: ");
    scanf("%lf", &x);

    double result = evaluateCubicSpline(x, data, coefficients, n);
    printf("Interpolated value at x = %lf is %lf\n", x, result);

    free(data);
    free(coefficients);

    return 0;
}

