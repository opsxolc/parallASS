#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"

#include <mpi.h>

#define PI 3.141592653589793238463
#define S 250 // chunk size
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
void master(double eps, int n, ChunkRand* cr, std::ofstream& error_file){
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

    error_file << "Init OK\n";
    std::flush(error_file);

    for (int i = 0; i < n; ++i){
        cr->get_random_chunk(chunks[i], S);
/*
        for (int j = 0; j < S; ++j)
            error_file << chunks[i][3*j] << " " << chunks[i][3*j + 1] << " " << chunks[i][3*j+2] << std::endl;
        std::flush(error_file);
*/
        MPI_Isend(chunks[i], S * 3, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &sends[i]);
        MPI_Irecv(&sums[i], 1, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD, &recvs[i]);
    }

    while (true) {
        int index;
        MPI_Testany(n, recvs, &index, &flag, &status);
        error_file << "TEST> " << index << ", " << flag << std::endl;
        if (flag) {
            sum += sums[index];
            p += S; // chunk size
            error_file << "SUCCESS: sol=" << P * sums[index] / S << ", sum=" << sums[index] << std::endl;
            error_file << "full: sum=" << sum << ", p=" << p << ", sol=" << P * sum / p << ", var=" << std::abs(SOLUTION - P * sum / p) << std::endl;
            if (std::abs(SOLUTION - P * sum / p) <= eps)
                break;

            if (!buf.empty()){
                error_file << "pop buf" << std::endl;
                copy_chunk(buf.pop(), chunks[index]);
            } else {
                error_file << "gen by myself" << std::endl;
                cr->get_random_chunk(chunks[index], S);
            }

            MPI_Isend(chunks[index], S * 3, MPI_DOUBLE, index + 1, 0, MPI_COMM_WORLD, &sends[index]);
            MPI_Irecv(&sums[index], 1, MPI_DOUBLE, index + 1, 1, MPI_COMM_WORLD, &recvs[index]);

            std::flush(error_file);
            break;
        }
/*
        MPI_Testany(n, sends, &index, &flag, &status);
        if (flag) {
            error_file << "send more\n";
            if (!buf.empty()){
                error_file << "pop buf" << std::endl;
                copy_chunk(buf.pop(), chunks[index]);
            } else {
                error_file << "gen by myself" << std::endl;
                cr->get_random_chunk(chunks[index], S);
            }

            MPI_Isend(chunks[index], S * 3, MPI_DOUBLE, index + 1, 0, MPI_COMM_WORLD, &sends[index]);
            break;
        }
*/

        error_file << "fill buf: " << buf.get_head()+1 << std::endl;
        buf.gen();
    }

    for (int i = 0; i < n; ++i) {
        free(chunks[i]);
    }

    free(chunks);
    free(sums);
    free(sends);
    free(recvs);

    error_file << "RESULT\nSolution=" << P * sum / p << ", Delta=" << std::abs(SOLUTION - P * sum / p) << ", Points=" << p << std::endl;
    std::flush(error_file);

    MPI_Abort(MPI_COMM_WORLD, 0);
}

void slave(int rank, std::ofstream& error_file){
    bool first_iteration = true;
    int flag;
    MPI_Status status;
    double* chunk = get_empty_chunk();
    double* sum_p = (double*) malloc(2 * sizeof(double));
    double sum;
    MPI_Request send, recv;

    while (true) {
        if (first_iteration)
            MPI_Irecv(chunk, S * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv);

        MPI_Wait(&recv, &status);
        if (status.MPI_ERROR != MPI_SUCCESS && status.MPI_ERROR != MPI_ERR_LASTCODE) {
            error_file << "Non success recv: " << status.MPI_ERROR << " count=" << status.count << std::endl;
            //return;
        } else {
            // error_file << "SUCCESS RECV!\n";
        }
        std::flush(error_file);

/*
        for (int j = 0; j < S; ++j)
            error_file << chunk[3*j] << " " << chunk[3*j + 1] << " " << chunk[3*j+2] << std::endl;
        std::flush(error_file);
*/

        MPI_Irecv(chunk, S * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv);

        //compute
        for (int j = 0; j < 1000; ++j) {
        sum = 0;
        for (int i = 0; i < S; ++i) {
            if (chunk[3*i+1] * chunk[3*i+1] + chunk[3*i+2] * chunk[3*i+2] <= 1) {
                sum += sqrt(chunk[3*i+1] * chunk[3*i+1] + chunk[3*i+2] * chunk[3*i+2]);
            }
        }
        }
        error_file << "sum=" << sum << std::endl;

        if (!first_iteration) {
            MPI_Wait(&send, &status);
            if (status.MPI_ERROR != MPI_SUCCESS && status.MPI_ERROR != MPI_ERR_LASTCODE) {
                error_file << "Non success send: " << status.MPI_ERROR << " count=" << status.count << std::endl;
                //return;
            } else {
                // error_file << "SUCCESS SEND!\n";
            }
            std::flush(error_file);
        }

        MPI_Isend(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send);

        if (first_iteration)
            first_iteration = false;
    }

    free(chunk);
}

int main(int argc, char **argv){
    // sqrt(x * x + z * z)
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
        // printf("Invalid number of processes %d, must be at least 2\n");
        std::cout << "Invalid number of processes " << size << ", must be at least 2" << std::endl;
        std::flush(std::cout);
        MPI_Finalize();
        return 0;
    }

    char file_name[30];
    sprintf(file_name, "error_%d.txt", rank);
    std::ofstream error_file;
    error_file.open(file_name);
    error_file << "Writing this to a file.\n";
    std::flush(error_file);

    // if (rank)
        // std::vector<double> a(n);

    // std::cout << "rank: " << rank << std::endl;

    ChunkRand cr(0, 2, -1, 1, -1, 1);

    if (!rank)
        master(eps, size - 1, &cr, error_file);
    else
        slave(rank, error_file);

    error_file.close();

    MPI_Finalize();

    return 0;
}
