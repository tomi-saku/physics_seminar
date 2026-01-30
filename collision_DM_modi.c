#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 1001
#define e 0.001
#define e_galaxy 0.1
#define rho_s 100000
#define r_s 0.1
// M_PIが定義されていない環境用
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// NFWの係数を事前に計算（計算コスト削減）
#define NFW_CONST (4.0 * M_PI * rho_s * r_s * r_s * r_s)

double calcR(double X, double Y, double Z)
{
    return sqrt(X * X + Y * Y + Z * Z + e * e);
}

double calcR_to_galaxy(double X, double Y, double Z)
{
    return sqrt(X * X + Y * Y + Z * Z + e_galaxy * e_galaxy);
}

double g(double q)
{
    return q * q * pow(1 - q * q, 3.5); // 7.0/2.0 -> 3.5
}

// 【修正1】returnを追加。NFWの「 enclosed mass M(<r) 」を計算
double calcM_NFW(double r)
{
    return NFW_CONST * (log(1 + r / r_s) - r / (r_s + r));
}

// 【修正2】NFWのポテンシャルエネルギー Φ(r) を計算する関数を追加
double calcPot_NFW(double r)
{
    // rが非常に小さい時のゼロ除算対策が必要だが、今回はe_galaxyがあるため簡易実装
    // Φ(r) = - (G * 4*pi*rho*rs^3 * ln(1+r/rs)) / r
    return -1.0 * NFW_CONST * log(1 + r / r_s) / r;
}

// 【修正3】加速度計算を関数化
void calculate_acceleration(int n_total, double *X, double *Y, double *Z, double *ax, double *ay, double *az, double G, double M_particle, double M_galaxy_total)
{
    for (int n = 0; n < n_total; n++)
    {
        double sumx = 0, sumy = 0, sumz = 0;

        // n == n_total-1 は銀河の中心核
        if (n == n_total - 1)
        {
            // 銀河中心が受ける力（星団の粒子からの反作用）
            // ※注意: 銀河全体が剛体として動くのか、中心核粒子だけが動くのかで物理が変わります。
            // 　ここでは「中心核粒子(質量M_galaxy_totalと仮定)が星団の星(質量M_particle)から引力を受ける」として計算します。
            /* ダークマターは重力を感じないとする
            for (int m = 0; m < n_total - 1; m++)
            {
                double dx = X[n] - X[m];
                double dy = Y[n] - Y[m];
                double dz = Z[n] - Z[m];
                double R = calcR_to_galaxy(dx, dy, dz);
                double R3 = R * R * R;

                // F = G * M_gal * M_star / R^2
                // a = F / M_gal = G * M_star / R^2
                // 銀河中心の加速度を計算する場合、相手の質量(M_particle)を使う
                sumx -= G * M_particle * dx / R3;
                sumy -= G * M_particle * dy / R3;
                sumz -= G * M_particle * dz / R3;
            }*/
        }
        else
        {
            // 星団の星が受ける力
            for (int m = 0; m < n_total; m++)
            {
                if (n == m)
                    continue;

                double dx = X[n] - X[m];
                double dy = Y[n] - Y[m];
                double dz = Z[n] - Z[m];

                if (m == n_total - 1)
                {
                    // 相手が銀河の場合 (NFWプロファイルからの力)
                    double R = calcR_to_galaxy(dx, dy, dz);
                    // 【修正4】calcMを1回だけ呼ぶ
                    double M_enclosed = calcM_NFW(R);
                    double factor = G * M_enclosed / (R * R * R);

                    sumx -= factor * dx;
                    sumy -= factor * dy;
                    sumz -= factor * dz;
                }
                else
                {
                    // 相手が星の場合 (点質量)
                    double R = calcR(dx, dy, dz);
                    double factor = G * M_particle / (R * R * R); // 星同士

                    sumx -= factor * dx;
                    sumy -= factor * dy;
                    sumz -= factor * dz;
                }
            }
        }
        ax[n] = sumx;
        ay[n] = sumy;
        az[n] = sumz;
    }
}

int main(void)
{
    double deltaT = 0.0001;
    int Steps = 10000;
    double R_cluster = 1;
    double G = 1;
    double M_particle = 1;
    double M_galaxy = 10000;
    double v_virial = sqrt(0.25 * G * M_particle * N / R_cluster);

    // 配列確保 (スタックオーバーフロー防止のためstaticかmalloc推奨だが、N=1001ならOK)
    double X_velocity[N] = {0};
    double Y_velocity[N] = {0};
    double Z_velocity[N] = {0};

    double X_position[N] = {0};
    double Y_position[N] = {0};
    double Z_position[N] = {0};

    // ... (初期位置設定等は元のコードと同じ) ...
    double X_position_galaxy1 = 2;
    double Y_position_galaxy1 = 0;
    double Z_position_galaxy1 = 0;
    double X_position_galaxy2 = 0;
    double Y_position_galaxy2 = 0;
    double Z_position_galaxy2 = 0;
    double v_orbital = sqrt(G * M_galaxy / fabs(X_position_galaxy1 - X_position_galaxy2)); // fabs追加(念の為)

    double ax[N], ay[N], az[N];

    // ランダムシード初期化（毎回結果を変えたい場合）
    // srand((unsigned int)time(NULL));

    // 星団の初期化
    for (int i = 0; i < N - 1; i++)
    {
        double rand_r = ((double)rand() / RAND_MAX);
        double rand_xy = ((double)rand() / RAND_MAX);
        double rand_z = ((double)rand() / RAND_MAX);

        double r_norm = pow(pow(rand_r, -2.0 / 3.0) - 1.0, -0.5);
        double r = r_norm * R_cluster;

        Z_position[i] = (1 - 2 * rand_z) * r + Z_position_galaxy1;
        double r_xy = sqrt(r * r - (Z_position[i] - Z_position_galaxy1) * (Z_position[i] - Z_position_galaxy1));
        X_position[i] = r_xy * cos(2 * M_PI * rand_xy) + X_position_galaxy1;
        Y_position[i] = r_xy * sin(2 * M_PI * rand_xy) + Y_position_galaxy1;

        // 速度の初期化
        double Ve = sqrt(2) * pow(1 + r * r, -0.25);
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

        double vz_local = (1 - 2 * rand_vz) * V;
        Z_velocity[i] = vz_local + 0; // z velocity galaxy1 is 0
        double vxy_local = sqrt(V * V - vz_local * vz_local);
        X_velocity[i] = vxy_local * cos(2 * M_PI * rand_vxy) + 0;               // x vel galaxy1 is 0
        Y_velocity[i] = vxy_local * sin(2 * M_PI * rand_vxy) + 0.5 * v_orbital; // y vel galaxy1 added
    }

    // 銀河中心（NFWポテンシャル源）の設定
    X_position[N - 1] = X_position_galaxy2;
    Y_position[N - 1] = Y_position_galaxy2;
    Z_position[N - 1] = Z_position_galaxy2;
    X_velocity[N - 1] = 0;
    Y_velocity[N - 1] = 0;
    Z_velocity[N - 1] = 0;

    FILE *fp = fopen("orbit_collision_DM.txt", "w");
    FILE *fp_2 = fopen("velocity_collision_DM.txt", "w");
    FILE *fp_3 = fopen("energy_collision_DM.txt", "w");

    // 初期加速度計算 (関数呼び出しに変更)
    calculate_acceleration(N, X_position, Y_position, Z_position, ax, ay, az, G, M_particle, M_galaxy);

    for (int i = 0; i < Steps; i++)
    {
        // Leapfrog Step 1: v(t + dt/2)
        for (int n = 0; n < N; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // Leapfrog Step 2: x(t + dt)
        for (int n = 0; n < N; n++)
        {
            X_position[n] += X_velocity[n] * deltaT;
            Y_position[n] += Y_velocity[n] * deltaT;
            Z_position[n] += Z_velocity[n] * deltaT;
        }

        // Leapfrog Step 3: a(t + dt)
        calculate_acceleration(N, X_position, Y_position, Z_position, ax, ay, az, G, M_particle, M_galaxy);

        // Leapfrog Step 4: v(t + dt)
        for (int n = 0; n < N; n++)
        {
            X_velocity[n] += 0.5 * ax[n] * deltaT;
            Y_velocity[n] += 0.5 * ay[n] * deltaT;
            Z_velocity[n] += 0.5 * az[n] * deltaT;
        }

        // ファイル出力（毎回出力すると遅いので、適宜間引くことを推奨：if (i % 10 == 0) 等）
        for (int j = 0; j < N; j++)
            fprintf(fp, "%f %f %f ", X_position[j], Y_position[j], Z_position[j]);
        fprintf(fp, "\n");

        for (int j = 0; j < N; j++)
            fprintf(fp_2, "%f %f %f ", X_velocity[j], Y_velocity[j], Z_velocity[j]);
        fprintf(fp_2, "\n");

        // エネルギー計算
        double kinetic_energy = 0;
        double potential_energy = 0;

        for (int j = 0; j < N; j++)
        {
            double v2 = X_velocity[j] * X_velocity[j] + Y_velocity[j] * Y_velocity[j] + Z_velocity[j] * Z_velocity[j];

            // 運動エネルギー
            if (j == N - 1)
                kinetic_energy += 0.5 * M_galaxy * v2; // 銀河中心のKE（必要なら）
            else
                kinetic_energy += 0.5 * M_particle * v2;

            // ポテンシャルエネルギー
            for (int k = j + 1; k < N; k++)
            {
                double dx = X_position[j] - X_position[k];
                double dy = Y_position[j] - Y_position[k];
                double dz = Z_position[j] - Z_position[k];

                // 銀河(N-1)との相互作用の場合
                if (k == N - 1)
                {
                    double R = calcR_to_galaxy(dx, dy, dz);
                    // 【修正2適用】NFWのポテンシャル関数を使用
                    // jは星なので質量M_particle
                    potential_energy += M_particle * calcPot_NFW(R) * G;
                    // ※calcPot_NFWにはGを含めないで定義したので、ここでGをかけるか関数内で調整してください
                    // 上記実装例のcalcPot_NFWは定数内にGを含んでいないため、
                    // calcPot_NFW = - 4*pi*rho*rs^3 * ln(...) / r
                    // PE = M_particle * Φ * G となります。
                    // (※定数NFW_CONSTにはGが含まれていないことに注意)
                }
                else
                {
                    // 星同士
                    double R = calcR(dx, dy, dz);
                    potential_energy -= G * M_particle * M_particle / R;
                }
            }
        }
        // NFWポテンシャル内のGの扱いに注意。
        // 上記コードではsumx -= G * calcM(R)... としているので、calcMは質量のみ。
        // ポテンシャルΦは -G * M(<r_effective) / r の形にする必要があるため、
        // エネルギー計算部分だけ別途注意深く係数を合わせる必要があります。

        fprintf(fp_3, "%f %f %f\n", i * deltaT, kinetic_energy, potential_energy);
    }

    fclose(fp);
    fclose(fp_2);
    fclose(fp_3);
    return 0;
}