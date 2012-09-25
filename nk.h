typedef struct {
  int n;
  int k;

  int **neighbors;
  int **friends;
  int *num_friends;
  int **is_neighbor;
  int *num_subproblems;
  int **subproblems;
  int *num_subproblems2;
  int **subproblems2;
  
  double **f;

  double optimum;
  
} NK_instance;

void generate_nk(int n, int k, NK_instance *nk);
double solve_nk(NK_instance *nk);
void bb_nk(int *x, int current, double current_f, int **index, int *index_size, double *max_contrib, double max_remain, NK_instance *nk,double ***best);

int num_subproblems_nk(NK_instance *nk, int low, int high);
void sort_int_array(int *x, int n);
double local_search_nk(NK_instance *nk);
void prepare_for_solve_nk(NK_instance *nk, int **index, int *index_size, double *max_contrib, double ***best);
double run_local_search_nk(char *x, NK_instance *nk,long *num_flips);
void load_nk(char *fname, NK_instance *nk);
void free_nk(NK_instance *nk);
