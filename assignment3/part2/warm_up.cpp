
void warm_up() {
    double dummy = 1.0;
    #pragma omp target data map(tofrom: dummy)
    {}
}