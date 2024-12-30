#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define N 50
void jacobi(double A_matrix[][N], double b[], double initial[]);
void seidel(double A_matrix[][N], double b[], double initial[]);
void sor(double A_matrix[][N], double b[], double initial[], double w);
void matrix_mult(double A_matrix[][N], double x[]);
int main() 
{
  double A_matrix [N][N];
  double b[N];
  double initial[N];
  double rowSum[N];
  double w = 1.1;
  double calc =0;
  char i = 0;
  char j = 0;
  char k = 0;
  clock_t start, end;
  //srand(time(NULL));
  for(i = 0; i<N;i++)
  {
    b[i] = i+1;
    initial[i] = 1;
  }
  for(i = 0; i<10; i++)
  {
    for(j = 0; j<N; j++)
    {
      for(k = 0; k<N; k++)
      {
        A_matrix[j][k] = (double) rand() / RAND_MAX;
        //printf("A_matrix[%d][%d] = %lf \n",j,k, A_matrix[j][k]);
        //printf("%lf" , RAND_MAX);
      }
    }
    for(j = 0; j< N ; j++)
    {
      for(k =0; k<N ; k++)
      {
        A_matrix[j][j]+= A_matrix[j][k];
      }  
    }
    printf("\n-----------------------------------matrix.%d-----------------------------------------\n",i+1);
    start= clock();
    jacobi(A_matrix, b, initial);
    end = clock();
    printf("\nThe time taken in jacobi for matrix %d = %lf microseconds\n", i+1 , (double)(end-start));
    start = clock();  
    seidel(A_matrix, b, initial);
    printf("\nThe time taken in seidel for matrix %d = %lf microseconds\n", i+1 , (double)(end-start));
    end= clock();
    for(j=0; j<5;j++)
    {
      w = 1 + ((double) rand() / RAND_MAX);
      printf("\nsor_%d_w = %lf\n",i+1,j+1, w);
      start = clock();
      sor(A_matrix, b, initial,w);
      end = clock();
      printf("end = %lf, start = %lf\n",end,start);
      printf("\nThe time taken in sor for matrix %d with w =%lf = %lf microseconds\n", i+1 , w, (double)(end-start));
    }  
  }
  
  
  return 0;
}

void jacobi(double A_matrix[][N], double b[], double initial[])
{
  //printf("Debug %lf\n", A_matrix[0][0]);
  double x[N];
  double calc = 0;
  double L_of_x = 0;
  double L_of_initial = 1;
  double fabs_x = 0;
  int count = 1;
  char i = 0;
  char j = 0;
  bool flag =true;
  while(flag ==true) //& count < 3)
  {
    for(i = 0;i<N;i++) 
    {
      calc = 0;
      for(j = 0;j<N;j++)
      {
        calc += A_matrix[i][j]*initial[j];
      }
      x[i] = (b[i]-calc + (A_matrix[i][i]*initial[i])) / A_matrix[i][i];      
      //printf("x[%d] = %lf\n",i, x[i]);   
    }
    L_of_x = 0;
    for(i = 0;i<N;i++)
    {
      fabs_x = fabs(x[i]);
      if(L_of_x < fabs_x)
      {
        L_of_x = fabs_x;
      }
      initial[i] = x[i];
    }
    printf("Iteration.%d; L_of_x = %lf\t",count, L_of_x);
    printf("L_of_initial = %lf\t", L_of_initial);
    printf("Difference = %lf\n", fabs(L_of_x - L_of_initial));
    if(fabs(L_of_x - L_of_initial) < 0.003)
    {
      printf("Iteration has been stopped\n");
      printf("The solution of jacobi ;\n");
      for(i = 0; i<N; i++)
      {
        printf("x[%d] = %lf  ", i+1 , x[i]);
        if(i+1 % 11 ==0)
        {
          printf("\n");
        }
        if(i+1 % 11 ==0)
        {
          printf("\n");
        }
       }
      matrix_mult(A_matrix,x);
      //printf("Iteration.%d and the final = %lf\n",count,fabs(L_of_x - L_of_initial));
      flag =false;
    }
    L_of_initial = L_of_x;
    count++;
  }
}

void seidel (double A_matrix[][N], double b[], double initial[])
{
  double x[N];
  double calc = 0;
  double L_of_x = 0;
  double L_of_initial = 1;
  double fabs_x = 0;
  int count = 1;
  char i = 0;
  char j = 0;
  bool flag = true;
  while(flag == true)
  {
    for(i = 0;i<N;i++)
    {
      calc = 0;
      for(j = 0;j<N;j++)
      {
        calc += A_matrix[i][j]*initial[j];
      }
      x[i] = (b[i]-calc + (A_matrix[i][i]*initial[i])) / A_matrix[i][i];
      //printf("x[%d] = %lf\n", i, x[i]);
      //printf("x[%d] = %lf\n", i, x[i]);
      //printf("initial[%d]= %lf\n",i, initial[i]);
      initial[i] = x[i];    
      //printf("initial[%d]= %lf\n",i, initial[i]);
    }
    L_of_x = 0;
    for(i = 0;i<N;i++)
    {
      fabs_x = fabs(x[i]);
      //printf("fabs_x[%d] = %lf\n ", i , fabs_x);
      if(L_of_x < fabs_x)
      {
        L_of_x = fabs_x;
      }  
    }
    printf("Iteration.%d; L_of_x = %lf\t",count, L_of_x);
    printf("L_of_initial = %lf\t", L_of_initial);
    printf("Difference = %lf\n", fabs(L_of_x - L_of_initial));
    if(fabs(L_of_x - L_of_initial) < 0.003)
    {
      printf("Iteration has been stopped\n");
      printf("The solution of siedel;\n");
      for(i = 0; i<N; i++)
      {
        printf("x[%d] = %lf  ", i+1 , x[i]);
        if(i+1 % 11 ==0)
        {
          printf("\n");
        }
          if(i+1 % 11 ==0)
          {
            printf("\n");
          }      
      }
      matrix_mult(A_matrix, x);
      //printf("Iteration.%d and the final = %lf\n",count,fabs(L_of_x - L_of_initial));
      flag =false;
    }    
    L_of_initial = L_of_x;
    count++;
  }  
}


void sor (double A_matrix[][N], double b[], double initial[], double w)
{
  double x[N];
  double calc = 0;
  double L_of_x = 0;
  double L_of_initial = 1;
  double fabs_x = 0;
  int count = 1;
  char i = 0;
  char j = 0;
  bool flag = true;
  while(flag == true)
  {
    for(i = 0;i<N;i++)
    {
      calc =0;
      for(j = 0;j<N;j++)
      {
        calc += A_matrix[i][j]*initial[j];
      }
      x[i] = (w*(b[i]-calc + (A_matrix[i][i]*initial[i])) / A_matrix[i][i]) + (1-w)*initial[i];
      //printf("x[%d] = %lf\n", i, x[i]);
      //printf("x[%d] = %lf\n", i, x[i]);
      //printf("initial[%d]= %lf\n",i, initial[i]);
      initial[i] = x[i];    
      //printf("initial[%d]= %lf\n",i, initial[i]);
    }
    L_of_x = 0;
    for(i = 0;i<N;i++)
    {
      fabs_x = fabs(x[i]);
      //printf("fabs_x[%d] = %lf\n ", i , fabs_x);
      if(L_of_x < fabs_x)
      {
        L_of_x = fabs_x;
      }  
    }
    printf("Iteration.%d; L_of_x = %lf\t",count, L_of_x);
    printf("L_of_initial = %lf\t", L_of_initial);
    printf("Difference = %lf\n", fabs(L_of_x - L_of_initial));
    if(fabs(L_of_x - L_of_initial) < 0.003)
    {
      printf("Iteration has been stopped\n");
      printf("The solution of sor with w = %lf;\n", w);
      for(i = 0; i<N; i++)
      {
        printf("x[%d] = %lf  ", i+1 , x[i]);
        if(i+1 % 11 ==0)
        {
          printf("\n");
        }
        if(i+1 % 11 ==0)
        {
          printf("\n");
        }
        
      }
      matrix_mult(A_matrix, x);
      //printf("Iteration.%d and the final = %lf\n",count,fabs(L_of_x - L_of_initial));
      flag =false;
    }
    L_of_initial = L_of_x;
    count++;
  }  
}

void matrix_mult(double A_matrix[][N] , double x[])
{
  char i =0;
  char j=0;
  double sum =0;
  printf("\n");
  for(i = 0;i<N;i++)
  {
    sum = 0;
    for(j = 0;j<N;j++)
    {
      sum += A_matrix[i][j];
    }
    printf("Ax[%d] = %lf ", i+1, sum*x[i]);
    if(i+1 % 11 == 0)
    {
    printf("\n");
    }
  } 
  
}
