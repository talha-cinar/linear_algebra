#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#define N 6
void jacobi(double A_matrix[][N], double b[], double initial[]);
void seidel(double A_matrix[][N], double b[], double initial[]);
void sor(double A_matrix[][N], double b[], double initial[], double w);
int main(void)
{
  double A_matrix[N][N]= {
                      {20.0, 5.0, 1.0, -1.0, 1.0, 5.0 },
                      {-2.0, 8.0, 0.0, -1.0, 0.0, 1.0},
                      {-5.0, 3.0, -21.0, -2.0, -4.0, 5.0},
                      {4.0, -2.0, 5.0, -25.0, -5.0, 0.0 },
                      {1.0, -4.0, 2.0, 1.0, 16.0, 4.0},
                      {2.0, -1.0, 5.0, 0.0,-4.0, -24.0}
  };

  double b[] = {-25, -68, 22, 12, 20, 67};
  double initial[] = {1,0,1,0,1,1};
  jacobi(A_matrix, b, initial);
  seidel(A_matrix, b, initial);
  sor(A_matrix, b, initial, 1.1);
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
      //printf("Iteration.%d and the final = %lf\n",count,fabs(L_of_x - L_of_initial));
      flag =false;
    }
    L_of_initial = L_of_x;
    count++;
  }  
}
