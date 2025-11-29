#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1000

double r=0;
double x=0;
double y=0;
double z=0;

double vx=0;
double vy=0;
double vz=0;

double R = 1;
double G = 1;
double M = 1;

double deltaR = 0.01;


int main(void){

double X_position[N] = {};
double Y_position[N] = {};
double Z_position[N] = {};

for(int i=0;i<N;i++){
double rand_r =((double)rand() / RAND_MAX);
double rand_xy =((double)rand() / RAND_MAX);
double rand_z =((double)rand() / RAND_MAX);
// double rand_vx =((double)rand() / RAND_MAX);
// double rand_vy =((double)rand() / RAND_MAX);
// double rand_vz =((double)rand() / RAND_MAX);

r = pow(pow(rand_r,-(double)2 / (double)3)-1,-(double)1/ (double)2);

Z_position[i]= (1-2*rand_z) * r;
X_position[i] = sqrt(r*r -z*z) * cos(2*M_PI*rand_xy);
Y_position[i] = sqrt(r*r -z*z) * sin(2*M_PI*rand_xy);
}
 FILE *fp = fopen("plummer_initial_position.txt", "w");
 for (int j = 0; j < N; j++) {
      fprintf(fp, "%f %f %f \n", X_position[j], Y_position[j], Z_position[j]);
    }
  fprintf(fp, "\n");
  fclose(fp);

  int count=0;
  int r_up=deltaR;
  for(int k; k<N;k++){
   
    int radius = sqrt(X_position[k]*X_position[k] + Y_position[k]*Y_position[k] + Z_position[k]*Z_position[k]);

    if(radius>deltaR)

  }

}