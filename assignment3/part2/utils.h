int coord_to_index(const int N, const double coord);

double index_to_coord(const int N, const int index);

double grid_spacing(const int N);

void init_edges(const int N, double ***a);

void subtract_arrays(const int N, double ***A, double ***B, double ***C);

int copy_grid_to_device(double* const host_pointer, double* const device_pointer, const int grid_size);

int copy_grid_from_device(double* const host_pointer, double* const device_pointer, const int grid_size);
