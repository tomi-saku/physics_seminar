#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1001
#define e 0.01
#define e_galaxy 0.1

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
    double G = 1;
    double M = 1;
    double M_galaxy = 100;
    double v_virial = sqrt(0.25 * G * M * N / R);

    double X_velocity[N] = {};
    double Y_velocity[N] = {};
    double Z_velocity[N] = {};

    double X_position[N] = {};
    double Y_position[N] = {};
    double Z_position[N] = {};

    double X_velocity_galaxy = 0;
    double Y_velocity_galaxy = 0;
    double Z_velocity_galaxy = 0;

    double X_position_galaxy = 2;
    double Y_position_galaxy = 0;
    double Z_position_galaxy = 0;

    // ✅ 加速度を保持する配列を追加
    double ax[N], ay[N], az[N];

    int plumer_r = 0;

    for (int i = 0; i < N - 1; i++)
    {
        while (true)
        {
            double rand_x = ((double)rand() / RAND_MAX);
            double rand_y = ((double)rand() / RAND_MAX);
            double rand_z = ((double)rand() / RAND_MAX);

            X_position[i] = 1 - 2 * rand_x;
            Y_position[i] = 1 - 2 * rand_y;
            Z_position[i] = 1 - 2 * rand_z;

            if (sqrt(X_position[i] * X_position[i] + Y_position[i] * Y_position[i] + Z_position[i] * Z_position[i]) <= R)
                break;
        }
        X_position[N - 1] = X_position_galaxy;
        Y_position[N - 1] = Y_position_galaxy;
        Z_position[N - 1] = Z_position_galaxy;

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
    X_position[N - 1] = X_position_galaxy;
    Y_position[N - 1] = Y_position_galaxy;
    Z_position[N - 1] = Z_position_galaxy;

    X_velocity[N - 1] = X_velocity_galaxy;
    Y_velocity[N - 1] = Y_velocity_galaxy;
    Z_velocity[N - 1] = Z_velocity_galaxy;

    FILE *fp = fopen("orbit_1000_collision.txt", "w");
    FILE *fp_2 = fopen("energy_1000_collision.txt", "w");

    // ✅ 初期加速度を計算
    for (int n = 0; n < N; n++)
    {
        double sumx = 0, sumy = 0, sumz = 0;
        if (n == N - 1)
        {
            for (int m = 0; m < N; m++)
            {
                if (m == n)
                    continue;
                double R = calcR_to_galaxy(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
                sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
                sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
            }
        }
        else
        {
            for (int m = 0; m < N; m++)
            {
                if (m == n)
                    continue;
                if (m == N - 1)
                {
                    double R = calcR_to_galaxy(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                    sumx -= G * M_galaxy * (X_position[n] - X_position[m]) / (R * R * R);
                    sumy -= G * M_galaxy * (Y_position[n] - Y_position[m]) / (R * R * R);
                    sumz -= G * M_galaxy * (Z_position[n] - Z_position[m]) / (R * R * R);
                }
                else
                {
                    double R = calcR(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                    sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
                    sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
                    sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
                }
            }
        }

        ax[n] = sumx;
        ay[n] = sumy;
        az[n] = sumz;
    }

    for (int i = 0; i < Steps; i++)
    {
        // ✅ Leapfrog: 速度を半ステップ更新
        for (int n = 0; n < N; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // ✅ 位置をフルステップ更新
        for (int n = 0; n < N; n++)
        {
            X_position[n] += X_velocity[n] * deltaT;
            Y_position[n] += Y_velocity[n] * deltaT;
            Z_position[n] += Z_velocity[n] * deltaT;
        }

        // ✅ 新しい加速度を再計算
        for (int n = 0; n < N; n++)
        {
            double sumx = 0, sumy = 0, sumz = 0;
            if (n == N - 1)
            {
                for (int m = 0; m < N; m++)
                {
                    if (m == n)
                        continue;
                    double R = calcR_to_galaxy(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                    sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
                    sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
                    sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
                }
            }
            else
            {
                for (int m = 0; m < N; m++)
                {
                    if (m == n)
                        continue;
                    if (m == N - 1)
                    {
                        double R = calcR_to_galaxy(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                        sumx -= G * M_galaxy * (X_position[n] - X_position[m]) / (R * R * R);
                        sumy -= G * M_galaxy * (Y_position[n] - Y_position[m]) / (R * R * R);
                        sumz -= G * M_galaxy * (Z_position[n] - Z_position[m]) / (R * R * R);
                    }
                    else
                    {
                        double R = calcR(X_position[n] - X_position[m], Y_position[n] - Y_position[m], Z_position[n] - Z_position[m]);
                        sumx -= G * M * (X_position[n] - X_position[m]) / (R * R * R);
                        sumy -= G * M * (Y_position[n] - Y_position[m]) / (R * R * R);
                        sumz -= G * M * (Z_position[n] - Z_position[m]) / (R * R * R);
                    }
                }
            }
            ax[n] = sumx;
            ay[n] = sumy;
            az[n] = sumz;
        }

        // ✅ 速度をもう半ステップ更新
        for (int n = 0; n < N; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // ファイル出力
        for (int j = 0; j < N; j++)
        {
            fprintf(fp, "%f %f %f ", X_position[j], Y_position[j], Z_position[j]);
        }
        fprintf(fp, "\n");

        double kinetic_energy = 0;
        double potential_energy = 0;

        for (int j = 0; j < N; j++)
        {
            if (j == N - 1)
            {
                kinetic_energy += M_galaxy * (X_velocity[j] * X_velocity[j] + Y_velocity[j] * Y_velocity[j] + Z_velocity[j] * Z_velocity[j]) / 2;
            }
            else
            {
                kinetic_energy += M * (X_velocity[j] * X_velocity[j] + Y_velocity[j] * Y_velocity[j] + Z_velocity[j] * Z_velocity[j]) / 2;
            }
            for (int k = j + 1; k < N; k++)
            {
                if (k == N - 1)
                {
                    double R = calcR_to_galaxy(X_position[j] - X_position[k], Y_position[j] - Y_position[k], Z_position[j] - Z_position[k]);
                    potential_energy += -G * M * M_galaxy / R;
                }
                else
                {
                    double R = calcR(X_position[j] - X_position[k], Y_position[j] - Y_position[k], Z_position[j] - Z_position[k]);
                    potential_energy += -G * M * M / R;
                }
            }
        }

        fprintf(fp_2, "%f %f %f\n", i * deltaT, kinetic_energy, potential_energy);
    }

    fclose(fp);
    fclose(fp_2);
    return 0;
}
