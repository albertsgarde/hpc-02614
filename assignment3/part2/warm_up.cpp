#include <omp.h>

void warm_up() {
    double dummy = 1.0;
    int dev = omp_get_default_device();
    #pragma omp target data map(tofrom: dummy) device(dev)
    {}
}