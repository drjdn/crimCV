#include <R.h>

void F77_SUB(srng)(void) {
  GetRNGstate();
}

void F77_SUB(erng)(void) {
  PutRNGstate();
}

double F77_SUB(runi)(void) {
  return unif_rand();
}
