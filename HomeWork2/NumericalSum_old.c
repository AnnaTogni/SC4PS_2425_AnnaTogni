//HomeWork2 -- 06-03-2026

// Addition of two vectors : d = a*x + y


// ******************************************************
// HOW TO COMPILE:

// gcc NumericalSum.C -o numericalsum
// ./numericalsum 5 2.0 
// in this case N = 5 and a=2.0

// otherwise 
// ./numericalsum
// and it will ask to insert N and a
// ******************************************************

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int N;
    double a;

    if (argc >= 3) {
        N = atoi(argv[1]);
        a = atof(argv[2]);
    } else {
        printf("Enter vector size N: ");
        scanf("%d", &N);

        printf("Enter constant a: ");
        scanf("%lf", &a);
    }

    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));

    char answer;
    double value;

    printf("Do you want all elements of x to be the same? (y/n): ");
    scanf(" %c", &answer);

    if (answer == 'y' || answer == 'Y') {
        printf("Value for x: ");
        scanf("%lf", &value);

        for (int i = 0; i < N; i++)
            x[i] = value;
    } else {
        for (int i = 0; i < N; i++) {
            printf("x[%d] = ", i);
            scanf("%lf", &x[i]);
        }
    }

    printf("Do you want all elements of y to be the same? (y/n): ");
    scanf(" %c", &answer);

    if (answer == 'y' || answer == 'Y') {
        printf("Value for y: ");
        scanf("%lf", &value);

        for (int i = 0; i < N; i++)
            y[i] = value;
    } else {
        for (int i = 0; i < N; i++) {
            printf("y[%d] = ", i);
            scanf("%lf", &y[i]);
        }
    }

    for (int i = 0; i < N; i++)
        d[i] = a * x[i] + y[i];

    printf("\nResult vector d:\n");
    for (int i = 0; i < N; i++)
        printf("d[%d] = %f\n", i, d[i]);

    free(x);
    free(y);
    free(d);

    return 0;
}