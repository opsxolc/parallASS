#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI 3.141592653589793238463
#define ITERATIONS 18

int nproc, rank;
double N, K;
double Lx, Ly, Lz;
double hx, hy, hz;
double tau;
double at;

int bi, ei, bj, ej, bk, ek;
int size_i, size_j, size_k;
int neighbors[6]; // i-1, i+1, j-1, j+1, k-1, k+1
MPI_Request sends[6], recvs[6];

class VecSUS {
    private: 
        double *v;

    public: 
        int size_i, size_j, size_k;

        VecSUS(){}

        VecSUS(int size_i, int size_j, int size_k) {
            this->size_i = size_i;
            this->size_j = size_j;
            this->size_k = size_k;

            v = (double *) malloc(size_i * size_j * size_k * sizeof(double));
        }
    
        double get(int i, int j, int k) {
            if (i * size_j * size_k + j * size_k + k < size_i * size_j * size_k)
                return v[i * size_j * size_k + j * size_k + k];

            std::cout << "oh shit we'r out of range in GET: " << i << " " << j << " " << k << std::endl;
            return -1;
        }

        void set(int i, int j, int k, double val) {
            if (i * size_j * size_k + j * size_k + k < size_i * size_j * size_k) {
                v[i * size_j * size_k + j * size_k + k] = val;
                return;
            }

            std::cout << "oh shit we'r out of range in SET\n"; 
        }

        double* operator[] (const int i) {
            if (0 < i && i < size_i) {
                return v + i * size_j * size_k;
            }
        }

        ~VecSUS() {
            free(v);
        }

};

VecSUS shadow_send_i, shadow_send_j, shadow_send_k;
VecSUS shadow_recv_i, shadow_recv_j, shadow_recv_k;

// TODO: local to global

double xi(int i) { 
    return hx * i;
}

double yj(int j) { 
    return hy * j;
}

double zk(int k) { 
    return hz * k;
}

double u(double x, double y, double z, double t) {
    return sin(PI/Lx * x)*sin(PI/Ly * y)*sin(PI/Lz * z)*cos(at * t);
}

double phi(double x, double y, double z) {
    return u(x, y, z, 0);
}

// depends on [i-1][j][k], [i][j-1][k], [i][j][k-1], [i+1][j][k], [i][j+1][k], [i][j][k+1]
double delta_h(VecSUS u, int i, int j, int k) {
    double dx = (u.get(i - 1, j, k) - 2 * u.get(i, j, k) + u.get(i + 1, j, k)) / (hx * hx); 
    double dy = (u.get(i, j - 1, k) - 2 * u.get(i, j, k) + u.get(i, j + 1, k)) / (hy * hy); 
    double dz = (u.get(i, j, k - 1) - 2 * u.get(i, j, k) + u.get(i, j, k + 1)) / (hz * hz); 

    return dx + dy + dz;
}

void set_send_shadows(VecSUS u){
    if (neighbors[0] > 0) {
        for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
                shadow_send_i.set(0, j, k, u.get(bi, j, k));
                if (neighbors[1] > 0)
                    shadow_send_i.set(1, j, k, u.get(ei, j, k));
            }
        }
    } if (neighbors[1] > 0) {
       for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
                shadow_send_i.set(0, j, k, u.get(ei, j, k));
            }
        } 
    }

     if (neighbors[2] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int k = bk; k <= ek; ++k) {
                shadow_send_j.set(0, i, k, u.get(i, bj, k));
                if (neighbors[3] > 0)
                    shadow_send_j.set(1, i, k, u.get(i, ej, k));
            }
        }
    } if (neighbors[3] > 0) {
       for (int i = bi; i <= ei; ++i) {
            for (int k = bk; k <= ek; ++k) {
                shadow_send_j.set(0, i, k, u.get(i, ej, k));
            }
        } 
    }

     if (neighbors[4] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int j = bj; j <= ej; ++j) {
                shadow_send_k.set(0, i, j, u.get(i, j, bk));
                if (neighbors[5] > 0)
                    shadow_send_k.set(1, i, j, u.get(i, j, ek));
            }
        }
    } if (neighbors[5] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int j = bj; j <= ej; ++j) {
                shadow_send_k.set(0, i, j, u.get(i, j, ek));
            }
        } 
    }
}

void set_recv_shadows(VecSUS u){
    if (neighbors[0] > 0) {
        for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
                u.set(bi, j, k, shadow_recv_i.get(0, j, k));
                if (neighbors[1] > 0)
                    u.set(ei, j, k, shadow_recv_i.get(1, j, k));
            }
        }
    } if (neighbors[1] > 0) {
       for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
                u.set(ei, j, k, shadow_recv_i.get(0, j, k));
            }
        } 
    }

     if (neighbors[2] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int k = bk; k <= ek; ++k) {
                u.set(i, bj, k, shadow_recv_j.get(0, i, k));
                if (neighbors[3] > 0)
                    u.set(i, ej, k, shadow_recv_j.get(1, i, k));
            }
        }
    } if (neighbors[3] > 0) {
       for (int i = bi; i <= ei; ++i) {
            for (int k = bk; k <= ek; ++k) {
                u.set(i, ej, k, shadow_recv_j.get(0, i, k));
            }
        } 
    }

     if (neighbors[4] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int j = bj; j <= ej; ++j) {
                u.set(i, j, bk, shadow_recv_k.get(0, i, j));
                if (neighbors[5] > 0)
                    u.set(i, j, ek, shadow_recv_k.get(1, i, j));
            }
        }
    } if (neighbors[5] > 0) {
        for (int i = bi; i <= ei; ++i) {
            for (int j = bj; j <= ej; ++j) {
                u.set(i, j, ek, shadow_recv_k.get(0, i, j));
            }
        } 
    }
}

void renew_shadow(VecSUS u) {
    set_send_shadows(u);

// TODO: add Irecv and wait for all

    if (neighbors[0] > 0) {
        MPI_Isend(shadow_send_i[0], shadow_send_i.size_j*shadow_send_i.size_k, 
            MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD, &sends[0]);
        if (neighbors[1] > 0) {
            MPI_Isend(shadow_send_i[1], shadow_send_i.size_j*shadow_send_i.size_k, 
                MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD, &sends[1]);
        }
    } else if (neighbors[1] > 0) {
        MPI_Isend(shadow_send_i[0], shadow_send_i.size_j*shadow_send_i.size_k, 
                MPI_DOUBLE, neighbors[1], 0, MPI_COMM_WORLD, &sends[1]);
    }

    if (neighbors[2] > 0) {
        MPI_Isend(shadow_send_j[0], shadow_send_j.size_j*shadow_send_j.size_k, 
            MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD, &sends[2]);
        if (neighbors[3] > 0) {
            MPI_Isend(shadow_send_j[1], shadow_send_j.size_j*shadow_send_j.size_k, 
                MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD, &sends[3]);
        }
    } else if (neighbors[3] > 0) {
        MPI_Isend(shadow_send_j[0], shadow_send_j.size_i*shadow_send_j.size_k, 
                MPI_DOUBLE, neighbors[3], 0, MPI_COMM_WORLD, &sends[3]);
    }

    if (neighbors[4] > 0) {
        MPI_Isend(shadow_send_k[0], shadow_send_k.size_j*shadow_send_k.size_k, 
            MPI_DOUBLE, neighbors[4], 0, MPI_COMM_WORLD, &sends[4]);
        if (neighbors[5] > 0) {
            MPI_Isend(shadow_send_k[1], shadow_send_k.size_j*shadow_send_k.size_k, 
                MPI_DOUBLE, neighbors[5], 0, MPI_COMM_WORLD, &sends[5]);
        }
    } else if (neighbors[5] > 0) {
        MPI_Isend(shadow_send_k[0], shadow_send_k.size_j*shadow_send_k.size_k, 
                MPI_DOUBLE, neighbors[5], 0, MPI_COMM_WORLD, &sends[5]);
    }
    
    set_recv_shadows(u);
}

void init(){
    hx = Lx / N;
    hy = Ly / N;
    hz = Lz / N;

    tau = 1 / K;

    at = PI * sqrt(1/(Lx*Lx) + 1/(Ly*Ly) + 1/(Lz*Lz));

    printf("hx=%lf hy=%lf hz=%lf tau=%lf at=%lf\nK=%lf\n", hx, hy, hz, tau, at, K);

    // TODO: init neighbors

    int shadow_i = (neighbors[0] > 0 ? 1 : 0) + (neighbors[1] > 0 ? 1 : 0);
    int shadow_j = (neighbors[2] > 0 ? 1 : 0) + (neighbors[3] > 0 ? 1 : 0);
    int shadow_k = (neighbors[4] > 0 ? 1 : 0) + (neighbors[5] > 0 ? 1 : 0);

    size_i = ei - bi + 1 + shadow_i;
    size_j = ej - bj + 1 + shadow_j;
    size_k = ek - bk + 1 + shadow_k;

    shadow_send_i = VecSUS(2, size_j-shadow_j, size_k-shadow_k);
    shadow_send_j = VecSUS(2, size_i-shadow_i, size_k-shadow_k);
    shadow_send_k = VecSUS(2, size_i-shadow_i, size_j-shadow_j);

    shadow_recv_i = VecSUS(2, size_j-shadow_j, size_k-shadow_k);
    shadow_recv_j = VecSUS(2, size_i-shadow_i, size_k-shadow_k);
    shadow_recv_k = VecSUS(2, size_i-shadow_i, size_j-shadow_j);
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

    init();

// u0
    VecSUS u0(size_i, size_j, size_k), u1(size_i, size_j, size_k);
    for (int i = bi; i <= ei; ++i) {
        for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
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

// TODO: renew shadow

// u1
    for (int i = bi; i <= ei; ++i) {
        for (int j = bj; j <= ej; ++j) {
            for (int k = bk; k <= ek; ++k) {
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
    for (int n = 0; n < ITERATIONS; ++n) { // 18 more steps
        diff = 0;
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                for (int k = 0; k <= N; ++k) {
                    if (i == 0 || j == 0 || k == 0 || i == N || j == N || k == N) {
                        vset(un_plus_1, i, j, k, 0);
                        continue;
                    }

                    // renew shadow
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
