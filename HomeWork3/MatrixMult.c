// Loop ordering for optimization

// Homework3 13-03-2026
// Computing the moltiplication between matrices (C = A x B)


// ******************************************************
// HOW TO COMPILE:

// gcc -O3 MatrixMult.c -o matrixmult
// gcc -O3 MatrixMult.c -lm -o matrixmult
// ./matrixmult output.txt

// ******************************************************


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// ******************************
// **** Function read matrix ****

void fill_matrix(double *A, int N, const char *name) {

    char answer;
    double value;

    // --- Case 1: a * Identity ---
    printf("Do you want %s = %c * I? (y/n): ", name, name[0] + 32);
    scanf(" %c", &answer);

    if (answer == 'y' || answer == 'Y') {

        printf("Enter value %c: ", name[0] + 32);
        scanf("%lf", &value);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i*N + j] = (i == j) ? value : 0.0;
            }
        }

    } else {

        // --- Case 2: Diagonal ---
        printf("Do you want %s to be diagonal? (y/n): ", name);
        scanf(" %c", &answer);

        if (answer == 'y' || answer == 'Y') {

            // initialize all to zero
            for (int i = 0; i < N*N; i++)
                A[i] = 0.0;

            printf("Enter %d diagonal elements of %s:\n", N, name);

            for (int i = 0; i < N; i++) {
                scanf("%lf", &A[i*N + i]);
            }

        } else {

            // --- Case 3: Full matrix ---
            printf("Enter elements of %s row by row:\n", name);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    scanf("%lf", &A[i*N + j]);
                }
            }
        }
    }
}

// ******************************

// *******************************************
// **** Function 2 : print multiplication ****

void print_three(double *A, double *B, double *C, int N) {

    for (int i = 0; i < N; i++) {

        // A
        for (int j = 0; j < N; j++)
            printf("%6.1f ", A[i*N + j]);

        if (i == N/2) printf("  x  ");
        else printf("     ");
        
        // B
        for (int j = 0; j < N; j++)
        printf("%6.1f ", B[i*N + j]);

        if (i == N/2) printf("  =  ");
        else printf("     ");
        
        // C (use optimized result)
        for (int j = 0; j < N; j++)
            printf("%6.1f ", C[i*N + j]);
            
        printf("\n");
    }
}

// *******************************************



int main(int argc, char *argv[]) {

    if (argc < 2) {
        printf("Usage: %s output_file\n", argv[0]);
        return 1;
    }

    char *filename = argv[1];

    printf("\n");
    
    int N;
    printf("Enter matrix size N: ");
    scanf("%d", &N);
    
    double *A = malloc(N * N * sizeof(double));
    double *B = malloc(N * N * sizeof(double));

    // Allocate two result matrices
    double (*C1)[N] = calloc(N, sizeof(*C1));
    double *C2 = calloc(N * N, sizeof(double));
    double (*C3)[N] = calloc(N, sizeof(*C3));
    double *C4 = calloc(N * N, sizeof(double));

    if (!A || !B || !C1 || !C2 || !C3 || !C4) {
        printf("Memory allocation failed\n");
        return 1;
    }

    printf("\nFill matrix A:\n");
    fill_matrix(A, N, "A");
    
    printf("\nFill matrix B:\n");
    fill_matrix(B, N, "B");

    double (*A1)[N] = (double (*)[N]) A;
    double (*B1)[N] = (double (*)[N]) B;


    printf("\n");

    
    // Benchmark start
    clock_t start1 = clock();

    // Naive Matrix Multiplication (C = A * B)
    // A, B, C are N x N matrices stored row-major

    for (int i = 0; i < N; i++) {
            printf("\rProgress: %d%%", (i * 100) / N);
            fflush(stdout);
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C1[i][j] += A1[i][k] * B1[k][j];
            }
        }
    }
    printf("\rProgress: 100%%\n");

    // Benchmark end
    clock_t end1 = clock();
    double time_spent1 = (double)(end1 - start1) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds -- NON OPTIMIZED - MATRIX FORMAT\n", time_spent1);


    // Benchmark start
    clock_t start2 = clock();

    // Naive Matrix Multiplication (C = A * B)
    // A, B, C are N x N matrices stored row-major

    for (int i = 0; i < N; i++) {
            printf("\rProgress: %d%%", (i * 100) / N);
            fflush(stdout);
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C2[i*N + j] += A[i*N + k] * B[k*N + j];
            }
        }
    }
    printf("\rProgress: 100%%\n");

    // Benchmark end
    clock_t end2 = clock();
    double time_spent2 = (double)(end2 - start2) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds -- NON OPTIMIZED - 1D ARRAY FORMAT\n", time_spent2);


    // Benchmark start
    clock_t start3 = clock();

    // Naive Matrix Multiplication (C = A * B)
    // A, B, C are N x N matrices stored row-major

    for (int i = 0; i < N; i++) {
            printf("\rProgress: %d%%", (i * 100) / N);
            fflush(stdout);
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                C3[i][j] += A1[i][k] * B1[k][j];
            }
        }
    }
    printf("\rProgress: 100%%\n");

    // Benchmark end
    clock_t end3 = clock();
    double time_spent3 = (double)(end3 - start3) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds -- OPTIMIZED - MATRIX FORMAT\n", time_spent3);



    // Benchmark start
    clock_t start4 = clock();

    // Optimized Matrix Multiplication (C = A * B)
    for (int i = 0; i < N; i++) {
            printf("\rProgress: %d%%", (i * 100) / N);
            fflush(stdout);
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                C4[i*N + j] += A[i*N + k] * B[k*N + j];
            }
        }
    }
    printf("\rProgress: 100%%\n");

    // Benchmark end
    clock_t end4 = clock();
    double time_spent4 = (double)(end4 - start4) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds -- OPTIMIZED - 1D ARRAY FORMAT\n", time_spent4);

    printf("\n");

    double times[4] = {time_spent1, time_spent2, time_spent3, time_spent4};
    double min = times[0];
    double max = times[0];
    for (int i = 1; i < 4; i++) {
        if (times[i] < min)
            min = times[i];
        if (times[i] > max)
            max = times[i];
    }
    double delta_time = max - min;
    printf("Time difference: %f seconds\n", delta_time);

    printf("\n");

    int equal = 1;

    for (int i = 0; i < N*N; i++) {
        if (fabs(C2[i] - C4[i]) > 1e-9) {
            equal = 0;
            break;
        }
    }

    if (equal){
        printf("The calculation has been successful!\n");
        if (N<6) 
            print_three(A, B, C4, N);
        }
    else
        printf("ERROR in the caltulation!");
    
    printf("\n");


    // Quick check for unitary and diagonal matrices 
    // double expected = N * a * b;

    // if (C[0] == expected)
    //     printf("Check PASSED: C elements = %f\n", expected);
    // else
    //     printf("Check FAILED\n");

    
    // Save matrix to file
    FILE *f = fopen(filename, "w");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(f, "%f ", C4[i*N + j]);
        }
        fprintf(f, "\n");
    }

    printf("\n");

    fclose(f);


    free(A);
    free(B);
    free(C1);
    free(C2);
    free(C3);

    return 0;
}