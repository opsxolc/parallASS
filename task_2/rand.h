class ChunkRand {
    private:
    double x1, x2, y1, y2, z1, z2;
    
    public:
    ChunkRand(double x1, double x2, double y1, double y2, double z1, double z2);
    void get_random_chunk(double* chunk, int size);

};
