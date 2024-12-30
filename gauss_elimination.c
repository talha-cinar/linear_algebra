#include <stdio.h>
#include <stdlib.h>

// Forward Declarations
void gaussian_elimination(float **A, float *b, int m, int n);
void back_substitution(float **A, float *b, float *x, int m, int n);

int main()
{
    int m, n, i, j;
    float **A, *b, *x;

    // Input the matrix size and elements
    printf("Enter the number of rows in the matrix: ");
    scanf("%d", &m);
    printf("Enter the number of columns in the matrix: ");
    scanf("%d", &n);

    // Dynamically allocate memory for A, b, and x
    A = (float **) malloc(m * sizeof(float *));
    b = (float *) malloc(m * sizeof(float));
    x = (float *) malloc(n * sizeof(float));
    for (i = 0; i < m; i++) {
        A[i] = (float *) malloc(n * sizeof(float));
    }

    printf("Enter the elements of the matrix: \n");
    for (i = 0; i < m; i++) 
    {
        for (j = 0; j < n; j++) 
        {
            if (scanf("%f", &A[i][j]) != 1) 
            {
                printf("Error: Invalid input.\n");
                exit(1); 
            }
        }
    }

    // Input the vector b
    printf("Enter the elements of the vector b: \n");
    for (i = 0; i < m; i++) 
    {
        if (scanf("%f", &b[i]) != 1) 
        {
            printf("Error: Invalid input.\n");
            exit(1);
        }
    }

    // Solve the system of equations
    gaussian_elimination(A, b, m, n);
    back_substitution(A, b, x, m, n);

    // Print the solutions
    printf("The solutions are:\n");
    for (i = 0; i < n; i++) 
    {
        printf("x%d = %f\n", i+1, x[i]);
    }

    // Free dynamically allocated memory
    for (i = 0; i < m; i++) 
    {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);

    return 0;
}

// Pivoting function
void pivot(float **A, float *b, int k, int m, int n) 
{
    int max_row = k;
    float max_val = A[k][k];

    // Find row with largest pivot
    for (int i = k+1; i < m; i++) 
    {
        if (abs(A[i][k]) > max_val) 
        {
            max_row = i;
            max_val = abs(A[i][k]);
        }
    }

    // Swap rows if necessary
    if (max_row != k) 
    {
        for (int j=0; j<n; j++)
        {
            float temp = A[k][j];
            A[k][j] = A[max_row][j];
            A[max_row][j] = temp;
        }
        float temp = b[k];
        b[k] = b[max_row];
        b[max_row] = temp;
    }
}

// Gaussian Elimination function
void gaussian_elimination(float **A, float *b, int m, int n)
{
    int i, j, k;

    for (k = 0; k < n-1; k++) 
    {
        pivot(A, b, k, m, n);
        for (i = k+1; i < m; i++) 
        {
            float factor = A[i][k] / A[k][k];
            b[i] = b[i] - factor * b[k];
            for (j = k; j < n; j++) 
            {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
        }
    }
}

// Back Substitution function
void back_substitution(float **A, float *b, float *x, int m, int n)
{
    int i, j;
    x[n-1] = b[n-1] / A[n-1][n-1];
    for (i = n-2; i >= 0; i--) 
    {
        x[i] = b[i];
        for (j = i+1; j < n; j++) 
        {
            x[i] = x[i] - A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }
}
