#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1000
#define N_another 200
#define e 0.01
#define e_galaxy 0.1
#define N_sum 1200

int count = 0;

double calcR(double X, double Y, double Z)
{
    double r = sqrt(X * X + Y * Y + Z * Z + e * e);
    return r;
}

double calcR_to_galaxy(double X, double Y, double Z)
{
    double r = sqrt(X * X + Y * Y + Z * Z + e_galaxy * e_galaxy);
    return r;
}

int main(void)
{
    double deltaT = 0.0001;
    int Steps = 10000;
    double R = 1;
    double R_another = 0.5;
    double G = 1;
    double M = 1;
    double M_galaxy = 100;
    double v_virial = sqrt(0.25 * G * M * N / R);
    double v_virial_another = sqrt(0.25 * G * M * N_another / R);

    double X_velocity[N_sum] = {};
    double Y_velocity[N_sum] = {};
    double Z_velocity[N_sum] = {};

    double X_position[N_sum] = {};
    double Y_position[N_sum] = {};
    double Z_position[N_sum] = {};

    double X_velocity_galaxy = 0;
    double Y_velocity_galaxy = 0;
    double Z_velocity_galaxy = 0;

    double X_position_galaxy = 2;
    double Y_position_galaxy = 0;
    double Z_position_galaxy = 0;

    // ✅ 加速度を保持する配列を追加
    double ax[N_sum], ay[N_sum], az[N_sum];

    int plumer_r = 0;

    // １個目の銀河の初期位置と初期速度を計算
    for (int i = 0; i < N; i++)
    {
        double rand_x = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_y = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_z = 1 - 2 * (double)rand() / RAND_MAX;

        double rand_l = sqrt(rand_x * rand_x + rand_y * rand_y + rand_z * rand_z);

        rand_x = rand_x / rand_l;
        rand_y = rand_y / rand_l;
        rand_z = rand_z / rand_l;

        X_position[i] = R * rand_x;
        Y_position[i] = R * rand_y;
        Z_position[i] = R * rand_z;

        double rand_vx = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_vy = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_vz = 1 - 2 * (double)rand() / RAND_MAX;

        double rand_vl = sqrt(rand_vx * rand_vx + rand_vy * rand_vy + rand_vz * rand_vz);

        rand_vx = rand_vx / rand_vl;
        rand_vy = rand_vy / rand_vl;
        rand_vz = rand_vz / rand_vl;

        X_velocity[i] = v_virial * rand_vx;
        Y_velocity[i] = v_virial * rand_vy;
        Z_velocity[i] = v_virial * rand_vz;
    }
    // ２個目の銀河の初期位置と初期速度を計算
    for (int i = N; i < N_sum; i++)
    {

        double rand_x = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_y = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_z = 1 - 2 * (double)rand() / RAND_MAX;

        double rand_l = sqrt(rand_x * rand_x + rand_y * rand_y + rand_z * rand_z);

        rand_x = rand_x / rand_l;
        rand_y = rand_y / rand_l;
        rand_z = rand_z / rand_l;

        X_position[i] = R_another * rand_x + X_position_galaxy;
        Y_position[i] = R_another * rand_y + Y_position_galaxy;
        Z_position[i] = R_another * rand_z + Z_position_galaxy;

        double rand_vx = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_vy = 1 - 2 * (double)rand() / RAND_MAX;
        double rand_vz = 1 - 2 * (double)rand() / RAND_MAX;

        double rand_vl = sqrt(rand_vx * rand_vx + rand_vy * rand_vy + rand_vz * rand_vz);

        rand_vx = rand_vx / rand_vl;
        rand_vy = rand_vy / rand_vl;
        rand_vz = rand_vz / rand_vl;

        X_velocity[i] = v_virial_another * rand_vx + X_velocity_galaxy;
        Y_velocity[i] = v_virial_another * rand_vy + Y_velocity_galaxy;
        Z_velocity[i] = v_virial_another * rand_vz + Z_velocity_galaxy;
    }

    FILE *fp = fopen("orbit_1000_200_collision.txt", "w");
    FILE *fp_2 = fopen("energy_1000_200_collision.txt", "w");

    // ✅ 初期加速度を計算
    for (int n = 0; n < N_sum; n++)
    {
        double sumx = 0, sumy = 0, sumz = 0;
        for (int m = 0; m < N_sum; m++)
        {
            if (m == n)
                continue;

            double R = calcR(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
            sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
            sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
            sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
        }

        ax[n] = sumx;
        ay[n] = sumy;
        az[n] = sumz;
    }

    for (int i = 0; i < Steps; i++)
    {
        // ✅ Leapfrog: 速度を半ステップ更新
        for (int n = 0; n < N_sum; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // ✅ 位置をフルステップ更新
        for (int n = 0; n < N_sum; n++)
        {
            X_position[n] += X_velocity[n] * deltaT;
            Y_position[n] += Y_velocity[n] * deltaT;
            Z_position[n] += Z_velocity[n] * deltaT;
        }

        // ✅ 新しい加速度を再計算
        for (int n = 0; n < N_sum; n++)
        {
            double sumx = 0, sumy = 0, sumz = 0;
            for (int m = 0; m < N_sum; m++)
            {
                if (m == n)
                    continue;

                double R = calcR(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
                sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
                sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
            }
            ax[n] = sumx;
            ay[n] = sumy;
            az[n] = sumz;
        }

        // ✅ 速度をもう半ステップ更新
        for (int n = 0; n < N_sum; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // ファイル出力
        for (int j = 0; j < N_sum; j++)
        {
            fprintf(fp, "%f %f %f ", X_position[j], Y_position[j], Z_position[j]);
        }
        fprintf(fp, "\n");

        double kinetic_energy = 0;
        double potential_energy = 0;

        for (int j = 0; j < N_sum; j++)
        {
            kinetic_energy += M * (X_velocity[j] * X_velocity[j] + Y_velocity[j] * Y_velocity[j] + Z_velocity[j] * Z_velocity[j]) / 2;

            for (int k = j + 1; k < N_sum; k++)
            {
                double R = calcR(X_position[j] - X_position[k], Y_position[j] - Y_position[k], Z_position[j] - Z_position[k]);
                potential_energy += -G * M * M / R;
            }
        }

        fprintf(fp_2, "%f %f %f\n", i * deltaT, kinetic_energy, potential_energy);
    }

    fclose(fp);
    fclose(fp_2);
    return 0;
}
