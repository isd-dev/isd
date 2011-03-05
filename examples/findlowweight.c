#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "isd.h"
#include "cmdline.h"

isd_parameters param;
isd_stats stats;
FILE *fout;

void print_word_found(unsigned long *current_word) {
  int i;
  fprintf(fout, "(");
  for (i = 0; i < param.n; ++i)
    fprintf(fout, "%d", isd_get_bit(current_word, i));
  fprintf(fout, ")\n");
}

int main(int argc, char *argv[]) {
  struct gengetopt_args_info args_info;
  unsigned long **matrix = NULL;
  int error_code;
  int c;
  int row_limb_size;
  int i, j;
  fout = stdout;
  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

  param.n = args_info.length_arg;
  param.k = args_info.dimension_arg;
  param.w = args_info.weight_arg;
  param.seed = args_info.seed_arg;
  param.px = args_info.px_arg;
  param.py = args_info.py_arg;
  param.m = args_info.m_arg;
  param.l = args_info.l_arg;
  param.c = args_info.c_arg;
  param.r = args_info.r_arg;
  param.kx = args_info.kx_arg;
  param.ky = args_info.ky_arg;
  param.max_iteration = args_info.max_iter_arg;
  param.max_words = args_info.max_found_arg;

  if (param.max_iteration == 0)
    param.max_iteration = (1UL << 31) - 1;

  param.callback = &print_word_found;

  /* isd_fprintf_parameters(fout, &param); */

  if (param.n * param.k == 0) {
    fprintf(fout, "()\n");
    exit(EXIT_SUCCESS);
  }

  row_limb_size = isd_size_in_limbs(param.n);

  matrix = (unsigned long **)malloc(param.k * sizeof(unsigned long *));
  matrix[0] = (unsigned long *)calloc(row_limb_size * param.k,
                                      sizeof(unsigned long));

  for (i = 1; i < param.k; ++i)
    matrix[i] = matrix[i-1] + row_limb_size;

  for (i = 0; i < param.k; ++i)
    for (j = 0; j < param.n; ++j) {
      do {
        c = getchar();
        if (c == EOF) {
          fprintf(stderr, "Wrong input, not enough '0' or '1' in the stream\n");
          free(matrix[0]);
          free(matrix);
          exit(EXIT_FAILURE);
        }
      } while ((c != '0') && (c != '1'));
      if (c == '1')
        isd_set_bit(matrix[i], j);
    }

  error_code = isd_low_weight_word((const unsigned long **)matrix, &param, &stats);

  free(matrix[0]);
  free(matrix);

  fprintf(stderr, "n_found: %d\n", stats.n_found);

  exit(EXIT_SUCCESS);
}
