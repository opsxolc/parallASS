#include <cstdlib>
#include "rand.h"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

ChunkRand::ChunkRand(double x1, double x2, double y1, double y2, double z1, double z2) {
    this->x1 = x1;
    this->x2 = x2;
    this->y1 = y1;
    this->y2 = y2;
    this->z1 = z1;
    this->z2 = z2;
    std::srand(322); // !
};

void ChunkRand::get_random_chunk(double* chunk, int size){
    if (size <= 0)
        return;

    for (int i = 0; i < size; ++i) {
        chunk[3 * i] = fRand(x1, x2);
        chunk[3 * i + 1] = fRand(y1, y2);
        chunk[3 * i + 2] = fRand(z1, z2);
    }
}
