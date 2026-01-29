#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1000

double r_norm, r;
double R = 1;

double g(double q)
{
  return q * q * pow(1 - q * q, (double)7 / (double)2);
}

int main(void)
{

  double X_position[N] = {};
  double Y_position[N] = {};
  double Z_position[N] = {};

  double X_velocity[N] = {};
  double Y_velocity[N] = {};
  double Z_velocity[N] = {};

  for (int i = 0; i < N; i++)
  {
    double rand_r = ((double)rand() / RAND_MAX);
    double rand_xy = ((double)rand() / RAND_MAX);
    double rand_z = ((double)rand() / RAND_MAX);

    r_norm = pow(pow(rand_r, -(double)2 / (double)3) - 1, -(double)1 / (double)2);
    r = r_norm * R;

    Z_position[i] = (1 - 2 * rand_z) * r;
    X_position[i] = sqrt(r * r - Z_position[i] * Z_position[i]) * cos(2 * M_PI * rand_xy);
    Y_position[i] = sqrt(r * r - Z_position[i] * Z_position[i]) * sin(2 * M_PI * rand_xy);

    double Ve = sqrt(2) * pow(1 + r * r, -(double)1 / (double)4);
    double q = 0;
    while (true)
    {
      double rand_q = ((double)rand() / RAND_MAX);
      double rand_p = ((double)rand() / RAND_MAX);

      if (rand_p * 0.1 < g(rand_q))
      {
        q = rand_q;
        break;
      }
    }
    double V = Ve * q;

    double rand_vxy = ((double)rand() / RAND_MAX);
    double rand_vz = ((double)rand() / RAND_MAX);

    Z_velocity[i] = (1 - 2 * rand_vz) * V;
    X_velocity[i] = sqrt(V * V - Z_velocity[i] * Z_velocity[i]) * cos(2 * M_PI * rand_vxy);
    Y_velocity[i] = sqrt(V * V - Z_velocity[i] * Z_velocity[i]) * sin(2 * M_PI * rand_vxy);
  }
  FILE *fp = fopen("plummer_initial_position.txt", "w");
  for (int j = 0; j < N; j++)
  {
    fprintf(fp, "%f %f %f \n", X_position[j], Y_position[j], Z_position[j]);
  }
  fprintf(fp, "\n");
  fclose(fp);
}