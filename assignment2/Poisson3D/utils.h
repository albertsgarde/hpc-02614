int coord_to_index(const int N, const double coord);

double index_to_coord(const int N, const int index);

double grid_spacing(const int N);

void init_edges(const int N, double ***a);

void subtract_arrays(const int N, double ***A, double ***B, double ***C);
