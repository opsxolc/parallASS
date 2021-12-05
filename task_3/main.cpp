#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"

#include <mpi.h>

#define PI 3.141592653589793238463
#define S 100 // chunk size
#define SOLUTION 4 * PI / 3
#define P 2*2*2

double* get_empty_chunk(){
    double* result = (double*) malloc(S * 3 * sizeof(double));
    return result;
}

class ChunkBuf {
    private:

    double** buf;
    int head, size;
    ChunkRand* cr;

    public:

    ChunkBuf(int size, ChunkRand* cr) {
        this->size = size;
        this->cr = cr;
        head = -1;
        buf = (double**) malloc(size * sizeof(double*));
        for (int i = 0; i < size; ++i)
            buf[i] = get_empty_chunk();
    }

    int get_head(){
        return head;
    }

    bool empty(){
        return head < 0;
    }

    bool full(){
        return head == size - 1;
    }

    void gen(){
        if (full())
            return;
        cr->get_random_chunk(buf[++head], S);
    }

    double* pop(){
        if (empty())
            return NULL;
        return buf[head--];
    }

    ~ChunkBuf(){
        for (int i = 0; i < size; ++i)
            free(buf[i]);
        free(buf);
    }

};

void copy_chunk(double* c1, double* c2){
    for (int i = 0; i < S * 3; ++i)
        c2[i] = c1[i];
}

// n - number of slaves
void master(double eps, int n, ChunkRand* cr, std::ostream& error_file){
    double startTime = MPI_Wtime();

    int flag;
    MPI_Status status;
    double sum = 0, p = 0; // number of points computed
    double** chunks = (double**) malloc(n * sizeof(double*));
    double* sums = (double*) malloc(n * sizeof(double));
    MPI_Request* sends = (MPI_Request*) malloc(n * sizeof(MPI_Request));
    MPI_Request* recvs = (MPI_Request*) malloc(n * sizeof(MPI_Request));
    ChunkBuf buf(n, cr);

    for (int i = 0; i < n; ++i) {
        chunks[i] = get_empty_chunk();
    }

    for (int i = 0; i < n; ++i){
        cr->get_random_chunk(chunks[i], S);

        MPI_Isend(chunks[i], S * 3, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &sends[i]);
        MPI_Irecv(&sums[i], 1, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD, &recvs[i]);
    }

    while (true) {
        int index;
        MPI_Testany(n, recvs, &index, &flag, &status);
        if (flag) {
            sum += sums[index];
            p += S; // chunk size
            // error_file << "full: sum=" << sum << ", p=" << p << ", sol=" << P * sum / p << ", var=" << std::abs(SOLUTION - P * sum / p) << std::endl;
            // std::flush(error_file);
            if (std::abs(SOLUTION - P * sum / p) <= eps)
                break;

            MPI_Test(&sends[index], &flag, &status);
            if (flag) {
                // error_file << "send> " << index << std::endl;

                if (!buf.empty()){
                    // error_file << "pop buf" << std::endl;
                    copy_chunk(buf.pop(), chunks[index]);
                } else {
                    // error_file << "gen by myself" << std::endl;
                    cr->get_random_chunk(chunks[index], S);
                }

                MPI_Isend(chunks[index], S * 3, MPI_DOUBLE, index + 1, 0, MPI_COMM_WORLD, &sends[index]);
            }

            MPI_Irecv(&sums[index], 1, MPI_DOUBLE, index + 1, 1, MPI_COMM_WORLD, &recvs[index]);

            continue;
        }

        MPI_Testany(n, sends, &index, &flag, &status);
        if (flag) {
            // error_file << "send more> " << index << std::endl;
            if (!buf.empty()){
                // error_file << "pop buf" << std::endl;
                copy_chunk(buf.pop(), chunks[index]);
            } else {
                // error_file << "gen by myself" << std::endl;
                cr->get_random_chunk(chunks[index], S);
            }

            MPI_Isend(chunks[index], S * 3, MPI_DOUBLE, index + 1, 0, MPI_COMM_WORLD, &sends[index]);
            continue;
        }

        // error_file << "fill buf: " << buf.get_head()+1 << std::endl;
        buf.gen();
        // std::flush(error_file);
    }

    for (int i = 0; i < n; ++i) {
        free(chunks[i]);
    }

    free(chunks);
    free(sums);
    free(sends);
    free(recvs);

    error_file << "RESULT\nSolution=" << P * sum / p << ", Delta=" << std::abs(SOLUTION - P * sum / p) << ", Points=" << p << std::endl;
    error_file << "Time=" << MPI_Wtime() - startTime << std::endl;
    std::flush(error_file);

    MPI_Abort(MPI_COMM_WORLD, 0);
}

// void slave(int rank, std::ofstream& error_file){
void slave(int rank){
    bool first_iteration = true;
    int flag;
    MPI_Status status;
    double *chunk = get_empty_chunk(), *chunk_recv = get_empty_chunk(), *tmp;
    double *sum_p = (double*) malloc(2 * sizeof(double));
    double sum, sum_send;
    MPI_Request send, recv;

    while (true) {
        if (first_iteration)
            MPI_Irecv(chunk_recv, S * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv);

        MPI_Wait(&recv, &status);

        tmp = chunk;
        chunk = chunk_recv;
        chunk_recv = tmp;

        MPI_Irecv(chunk_recv, S * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv);

        //compute
        // for (int j = 0; j < 500; ++j) {
        sum = 0;
        for (int i = 0; i < S; ++i) {
            if (chunk[3*i+1] * chunk[3*i+1] + chunk[3*i+2] * chunk[3*i+2] <= 1) {
                sum += sqrt(chunk[3*i+1] * chunk[3*i+1] + chunk[3*i+2] * chunk[3*i+2]);
            }
        }
        // }

        if (!first_iteration) {
            MPI_Wait(&send, &status);
            if (status.MPI_ERROR != MPI_SUCCESS && status.MPI_ERROR != MPI_ERR_LASTCODE) {
            } else {
            }
        }

        sum_send = sum;

        MPI_Isend(&sum_send, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send);

        if (first_iteration)
            first_iteration = false;
    }

    free(chunk);
}

int main(int argc, char **argv){
    // sqrt(y * y + z * z)
    // 0 <= x <= 2, y * y + z * z <= 1

    if (argc != 2) {
        std::cout << "Invalid number of parameters " << argc << ", must be 2 (programm, eps)\n";
    }
    
    double eps;
    eps = strtod(argv[1], NULL);

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

    ChunkRand cr(0, 2, -1, 1, -1, 1);

    if (!rank)
        master(eps, size - 1, &cr, std::cout);
    else
        slave(rank);

    MPI_Finalize();

    return 0;
}
