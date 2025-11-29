
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1000
#define e 0.1

int count = 0;

  double calcR(double X,double Y,double Z){

  //イプシロンを導入する。これがないと、２粒子間の距離が小さくなった時に1/Rが発散する。
  double r =sqrt(X*X + Y*Y + Z*Z + e*e);

  return r;
  }

int main(void) {

  double deltaT = 0.001;
  int Steps =1000;
  double R = 1;
  double G =1;
  double M = 1;
  double v_virial = sqrt(G * M / R);

  double X_velocity[N]= {};
  double Y_velocity[N]= {};
  double Z_velocity[N]= {};

  double X_position[N]= {};
  double Y_position[N]= {};
  double Z_position[N]= {};

  double x = 1.0;
  double y = 1.0;
  double z = 1.0;

  double kinetic_energy = 0;
  double potential_energy = 0;

  for(int i = 0 ; i<N ; i++){

    while(true){
     double rand_x = ((double)rand()/RAND_MAX);
     double rand_y = ((double)rand()/RAND_MAX);
     double rand_z = ((double)rand()/RAND_MAX);

     // x = rand_x/1000;
     // y = rand_y/1000;
     // z = rand_z/1000;
     
     
    X_position[i] = 1 - 2 * rand_x;
    Y_position[i] = 1 - 2 * rand_y;
    Z_position[i] = 1 - 2 * rand_z;

     if(sqrt(X_position[i]*X_position[i] + Y_position[i]*Y_position[i] + Z_position[i]*Z_position[i]) <= R){
       break;
     }
    };

    while(true){
    double rand_vx = (double)rand()/RAND_MAX;
    double rand_vy = (double)rand()/RAND_MAX;
    double rand_vz = (double)rand()/RAND_MAX;

    X_velocity[i] = v_virial/10*(1 - 2*rand_vx);
    Y_velocity[i] = v_virial/10*(1 - 2*rand_vy);
    Z_velocity[i] = v_virial/10*(1 - 2*rand_vz);

    if(sqrt(X_velocity[i]*X_velocity[i] + Y_velocity[i]*Y_velocity[i] + Z_velocity[i]*Z_velocity[i]) <= v_virial/10){
       break;
     }

    }

    
  
  }

 FILE *fp = fopen("orbit_1000.txt","w");
 FILE *fp_2 = fopen("energy_1000.txt","w");

  for (int i = 0; i < Steps; i++) {

    for (int n = 0; n < N; n++) {

    double ax = 0;
    double ay = 0;
    double az = 0;

    for (int m = 0; m<N ; m++){
      if (m == n) continue; 
      double R = calcR(X_position[n]-X_position[m],Y_position[n]-Y_position[m],Z_position[n]-Z_position[m]);

      ax = ax -(G * M * (X_position[n]-X_position[m]) / (R*R*R));
      ay = ay -(G * M * (Y_position[n]-Y_position[m]) / (R*R*R));
      az = az -(G * M * (Z_position[n]-Z_position[m]) / (R*R*R));

    }


    X_velocity[n] += ax * deltaT;
    Y_velocity[n] += ay * deltaT;
    Z_velocity[n] += az * deltaT;

    X_position[n] += X_velocity[n] * deltaT;
    Y_position[n] += Y_velocity[n] * deltaT;
    Z_position[n] += Z_velocity[n] * deltaT;


    }

 //数値解をファイルに書き込む
    for(int j = 0 ; j<N ; j++){
      fprintf(fp,"%f %f %f " ,X_position[j],Y_position[j],Z_position[j]);
      // printf("%f %f %f " ,X_position[j],Y_position[j],Z_position[j]);
    }
      fprintf(fp,"\n");
      // printf("\n");

    double kinetic_energy = 0;
    double potential_energy = 0;
    
    for(int j = 0 ; j<N ; j++){
      count++;
      kinetic_energy += M * (X_velocity[j]*X_velocity[j] + Y_velocity[j]*Y_velocity[j] + Z_velocity[j]*Z_velocity[j]) /2;
      for(int k = j+1 ; k<N ; k++){
         double R = calcR(X_position[j]-X_position[k],Y_position[j]-Y_position[k],Z_position[j]-Z_position[k]);
         
      potential_energy += -G * M * M /R ;
      }
    }
      fprintf(fp_2,"%f %f %f ", i*deltaT,kinetic_energy,potential_energy);
      
      fprintf(fp_2,"\n");
printf("%d\n",count);

  };
  fclose(fp);
  fclose(fp_2);
  return 0;
};

