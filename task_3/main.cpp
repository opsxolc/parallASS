#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI 3.141592653589793238463

double N, K;
double Lx, Ly, Lz;
double hx, hy, hz;
double tau;
double at;

double xi(int i) { 
    return hx * i;
}

double yj(int j) { 
    return hy * j;
}

double zk(int k) { 
    return hz * k;
}

double vget(std::vector<double> v, int i, int j, int k) {
    if (i * (N+1) * (N+1) + j * (N+1) + k < (N+1) * (N+1) * (N+1))
        return v[i * (N+1) * (N+1) + j * (N+1) + k];

    std::cout << "oh shit we'r out of range in GET: " << i << " " << j << " " << k << std::endl;
    return -1;
}

void vset(std::vector<double> &v, int i, int j, int k, double val) {
    if (i * (N+1) * (N+1) + j * (N+1) + k < (N+1) * (N+1) * (N+1)) {
        v[i * (N+1) * (N+1) + j * (N+1) + k] = val;
        return;
    }

    std::cout << "oh shit we'r out of range in SET\n"; 
}

double u(double x, double y, double z, double t) {
    return sin(PI/Lx * x)*sin(PI/Ly * y)*sin(PI/Lz * z)*cos(at * t);
}

double phi(double x, double y, double z) {
    return u(x, y, z, 0);
}

double delta_h(std::vector<double> u, int i, int j, int k) {
    double dx = (vget(u, i - 1, j, k) - 2 * vget(u, i, j, k) + vget(u, i + 1, j, k)) / (hx * hx); 
    double dy = (vget(u, i, j - 1, k) - 2 * vget(u, i, j, k) + vget(u, i, j + 1, k)) / (hy * hy); 
    double dz = (vget(u, i, j, k - 1) - 2 * vget(u, i, j, k) + vget(u, i, j, k + 1)) / (hz * hz); 

    return dx + dy + dz;
}

int main(int argc, char **argv){

    if (argc != 6) {
        std::cout << "Invalid number of parameters " << argc << ", must be 6 (program, Lx, Ly, Lz, N, K)\n";
    }
    
    Lx = strtod(argv[1], NULL);
    Ly = strtod(argv[2], NULL);
    Lz = strtod(argv[3], NULL);

    N = strtod(argv[4], NULL);
    K = strtod(argv[5], NULL);

    hx = Lx / N;
    hy = Ly / N;
    hz = Lz / N;

    tau = 1 / K;

    at = PI * sqrt(1/(Lx*Lx) + 1/(Ly*Ly) + 1/(Lz*Lz));

    printf("hx=%lf hy=%lf hz=%lf tau=%lf at=%lf\nK=%lf\n", hx, hy, hz, tau, at, K);


// u0
    std::vector<double> u0((N+1)*(N+1)*(N+1)), u1((N+1)*(N+1)*(N+1));
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                if (i == 0 || j == 0 || k == 0 || i == N || j == N || k == N) {
                    vset(u0, i, j, k, 0);
                    continue;
                }

                double x_i = xi(i), y_j = yj(j), z_k = zk(k);
                
                double u0_ijk = phi(x_i, y_j, z_k);
                vset(u0, i, j, k, u0_ijk);
            }
        }
    }

// u1
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                if (i == 0 || j == 0 || k == 0 || i == N || j == N || k == N) {
                    vset(u1, i, j, k, 0);
                    continue;
                }

                vset(u1, i, j, k, vget(u0, i, j, k) + tau * tau / 2 * delta_h(u0, i, j, k));
            }
        }
    }

// un
    std::vector<double> un_plus_1((N+1)*(N+1)*(N+1)), un_sub_1 = u0, un = u1;
    double diff;
    for (int n = 0; n < 18; ++n) { // 18 more steps
        diff = 0;
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    if (i == 0 || j == 0 || k == 0 || i == N || j == N || k == N) {
                        vset(un_plus_1, i, j, k, 0);
                        continue;
                    }

                    double un_plus_1_ijk = tau * tau * delta_h(un, i, j, k) + 2 * vget(un, i, j, k) - vget(un_sub_1, i, j, k);

                    vset(un_plus_1, i, j, k, un_plus_1_ijk);

                    diff = fmax(diff, fabs(un_plus_1_ijk - u(xi(i), yj(j), zk(k), tau * (n + 2))));
                    // if (fabs(un_plus_1_ijk - u(xi(i), yj(j), zk(k), tau * (n + 2))) > 0)
                        // printf("(%d, %d, %d) = %f - %f\n", i, j, k, un_plus_1_ijk, u(xi(i), yj(j), zk(k), tau * (n + 2)));
                }
            }
        }
        un_sub_1 = un;
        un = un_plus_1;
    }

// output

    std::cout << "diff=" << diff << std::endl;

    // printf("\n\n");

    // for (int i = 0; i <= N; ++i) {
    //     for (int j = 0; j <= N; ++j) {
    //         printf("u1[%d][%d] = [", i, j);
    //         for (int k = 0; k <= N; ++k) {
    //             if (k != N)
    //                 printf("%.2f, ", vget(un, i, j, k));
    //             else
    //                 printf("%.2f", vget(un, i, j, k));
    //         }
    //         printf("]\n");
    //     }
    // }

    /*
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0 && size < 2) {
        std::cout << "Invalid number of processes " << size << ", must be at least 2" << std::endl;
        std::flush(std::cout);
        MPI_Finalize();
        return 0;
    }
    */

    // MPI_Finalize();

    return 0;
}
