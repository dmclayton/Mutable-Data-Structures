// bf-correct-attk.c
//
// Compile with gcc -std=c11 -O3.
//
// Simulate the attack on classic BFs described in Section 4.1. The attack has
// the following parameters:
//   k - number of hashes
//   m - 2^m = number of filter bits
//   n - size of the set
//   r - number of error queries
//   q - 2^q = number of test queries
//
// The adversary chooses a set of 2^q test queries and a set of r target queries.
// It searches for a subset of n test queries for which each of the target
// queries are a false positive of the Bloom filter representation.
//
// The inputs are simulated, meaning for each query, we choose k elements of
// [1..2^m] uniformly and independently, and run the attack using these values.
//
// The output of this program are the success rate and a measure of average-case
// time complexity of the attack for 1000 simulations, and varying values of the
// parameters.
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"

// Executes the attack for target queries encoded by r_filt and test queries
// encoded by t_filts for parameters k, m, n, r, and q, and returns 1 if the
// attack succeeded and 0 otherwise.
//
// Inputs:
//
//  r_filt - a Bloom filter of the target set.
//
//  t_filts - a Bloom filter for each element of the test set.
//
//  n_filts, st - state for carrying out the computation. This is here in order
//    to speed up the simulation.
//
//  k, m, n, r, q - attack parameters.
//
//  dbg - if this is > 0, then attack() outputs a trace of the computation for
//    debugging.
//
//  test_ct - used to keep track of the number of sets tested as the attack is
//    carried out.
int attack(char *r_filt, char **t_filts, char **n_filts, int *st,
           int k, int m, int n, int r, int q, int dbg, int *test_ct);

// Executes a simpler attack whereby the adversary tests randomly chosen subsets
// of n test queries until one is found that covers the target set. returns 1 if
// the attack succeeded and 0 otherwise.
//
// Inputs:
//  r_filt, t_filts - as above
//  k, m, n, r, q - as above
//  t - the maximum number of sampled subsets tried before giving up.
int ref_attack(char *r_filt, char **t_filts,
               int k, int m, int n, int r, int q, int t);

// Tests attack() against ref_attack() on simulated inputs. Returns 1 if
// both succeed, both fail, or attack() succeeds and ref_attack() fails.
//
// Inputs:
//  k, m, n, r, q - as above.
//  dbg - if this is > 0, then outputs the attack() trace.
//  seed - used to seed srand(). This is helpful for debugging.
//
// NOTE t=1000 is hard-coded.
//
// NOTE It's likely that attack() should succeed while ref_attack() fails, since
// the latter fails with some probability. However, if attack() fails and
// ref_attack() passes, then this is indicates a problem with attack().
int test(int k, int m, int n, int r, int q, int dbg, int seed);

// Executes test() a number of times and returns 1 if and only if all tests
// passed. If a test fails, then a message is output to stderr with the input
// seed.
//
// Inputs:
//  k, m, n, r, q - as above
//  trials - number of times Test() is run.
int multi_test(int k, int m, int n, int r, int q, int trials);

// Executes attack() a number of times. It sets avg_tests to be the average
// number of sets tested by attack(); that is, the average value of test_ct when
// the function returns. Returns the fraction of simulations in which the attack
// succeeded.
//
// Inputs:
//  k, m, n, r, q, trials - above
//  avg_tests - average value of test_ct over all trials
double sim(int k, int m, int n, int r, int q, int trials, double *avg_tests);

// Sets the idx-th bit of filter to 1.
void set_idx(int idx, char *filter) {
  int byte = idx >> 3;
  int bit = (byte << 3) ^ idx;
  filter[byte] |= 1 << bit;
}

// Sets the idx-th bit of filter to 0.
void unset_idx(int idx, char *filter) {
  int byte = idx >> 3;
  int bit = (byte << 3) ^ idx;
  filter[byte] &= ~(1 << bit);
}

// Returns the value of the idx-th bit of filter.
int get_idx(int idx, const char *filter) {
  int byte = idx >> 3;
  int bit = (byte << 3) ^ idx;
  return (filter[byte] & (1 << bit)) > 0;
}

static const int mask[] = {
  0x00000001, 0x00000003, 0x00000007, 0x0000000f,
  0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff,
  0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff,
  0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff,
  0x0001ffff, 0x0003ffff, 0x0007ffff, 0x000fffff,
  0x001fffff, 0x003fffff, 0x007fffff, 0x00ffffff,
  0x01ffffff, 0x03ffffff, 0x07ffffff, 0x0fffffff,
  0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};

// Returns a random index into filter of length 2^m.
int rand_idx(int m) {
  return rand() & mask[m-1];
}

// Generates a random, 2^m-bit Bloom filter with n elements and k hashes.
void random_filter(char *filter, int k, int m, int n) {
  int filter_bits = 1<<m;
  memset(filter, 0, filter_bits / 8 * sizeof(char));
  for (int i = 0; i < n*k; i++) {
    set_idx(rand_idx(m), filter);
  }
}

// Prints a 2^m-bit filter.
void print_filter(const char *filter, int m) {
  int filter_bits = 1<<m;
  for (int idx = 0; idx < filter_bits; idx++) {
    printf("%d", get_idx(idx, filter));
  }
  printf("\n");
}

// Swaps random elements of array idx of length 2^q.
void shuffle(int *idx, int q) {
  int j;
  for (int i = 0; i < 1<<q; i++) {
    j = rand_idx(q);
    idx[i] ^= idx[j];
    idx[j] ^= idx[i];
    idx[i] ^= idx[j];
  }
}

// Implements ref_attack().
int ref_attack(char *r_filt, char **t_filts,
               int k, int m, int n, int r, int q, int t) {
  int filter_bits = 1 << m;
  int filter_bytes = filter_bits / 8;
  int num_queries = 1 << q;

  // Initialize index. This is used for choosing random subsets of trial
  // queries.
  int *idx = malloc(num_queries * sizeof(int));
  for (int i = 0; i < num_queries; i++) {
    idx[i] = i;
  }

  char *m_filt = malloc(filter_bytes * sizeof(char));

  // Run attack.
  int res;
  for (int ct = 0; ct < t; ct++) {

    shuffle(idx, q);
    for (int j = 0; j < filter_bytes; j++) {
      m_filt[j] = ~r_filt[j];
      for (int i = 0; i < n; i++) {
        m_filt[j] |= t_filts[idx[i]][j];
      }
    }

    res = 1;
    for (int j = 0; j < filter_bytes; j++) {
      res = res && ((char)m_filt[j] == (char)0xff);
    }

    if (res) { // Success!
      break;
    }
  }
  free(m_filt);
  return res;
}

// Recursively searches for a solution, called by attack().
//
// Inputs:
//
//  r_filt, t_filts - as above
//
//  n_filts - a stack of Bloom filters representing the test sets in the path to
//    the root of the tree from the current node.
//
//  st - a stack of indices into t_filts corresponding to the test sets in the
//    path to the root.
//
//  filter_bytes - length of filter (in bytes, i.e. (2^m)/8)
//
//  num_queries - number of queries (2^q)
//
//  n - as above
//
//  h - the height of the stack (i.e. depth of the current node in the tree)
//
//  i - the last bit covered
//
//  dbg, test_ct - as above
int attack0(char *r_filt, char **t_filts, char **n_filts, int *st,
            int filter_bytes, int num_queries, int n, int h, int i, int dbg,
            int *test_ct) {
  // Max depth reached.
  if (h > n+1) {
    return 0;
  }
  if (dbg) {
    for (int x = 0; x < h; x++) {
      printf("%d ", st[x]);
    }
    printf("\n");
  }

  // Add 1 to the number of tested sets.
  if (test_ct != NULL) {
    (*test_ct)++;
  }

  // Check if r_filt = r_filt & n_filts[h]; if so, success!
  int res = 1;
  for (int ell = 0; ell < filter_bytes; ell++) {
    res = res && ((char)(~r_filt[ell] | n_filts[h-1][ell]) == (char)0xff);
  }
  if (res) {
    return 1;
  }

  // Find the next uncovered bit and look for a t_filt[t] that covers it that is
  // not on the stack.
  for (int j = i; j < (filter_bytes << 3); j++) {

    if (get_idx(j, r_filt) && !get_idx(j, n_filts[h-1])) {

      for (int t = 0; t < num_queries; t++) {

        // Check if t is on the stack.
        int ok = 1;
        for (int x = 1; x < h; x++) {
          ok = ok && (t != st[x]);
        }

        if (ok && get_idx(j, t_filts[t])) {

          // Compute new n_filts[h].
          for (int ell = 0; ell < filter_bytes; ell++) {
            n_filts[h][ell] = n_filts[h-1][ell] | t_filts[t][ell];
          }
          st[h] = t;
          res = attack0(r_filt, t_filts, n_filts, st, filter_bytes, num_queries,
                        n, h+1, j, dbg, test_ct);
          if (res) {
            return 1;
          }
        }
      }
    }
  }
  return 0;
}

int attack(char *r_filt, char **t_filts, char **n_filts, int *st,
           int k, int m, int n, int r, int q, int dbg, int *test_ct) {
  int filter_bits = 1 << m;
  int filter_bytes = filter_bits >> 3;
  int num_queries = 1 << q;

  // Quickly check if there's no solution.
  //
  // If this check passes, this does not imply there is a solution. This code is
  // just to eliminate a lot of worst-case searches.
  for (int j = 0; j < filter_bits; j++) {
    if (get_idx(j, r_filt)) {
      int res = 0;
      for (int t = 0; t < num_queries; t++) {
        res = res || (get_idx(j, t_filts[t]));
      }
      if (!res) {
        return 0;
      }
    }
  }

  // Execute attack.
  for (int i = 0; i < num_queries; i++) {
    memset(n_filts[i], 0, filter_bytes * sizeof(char));
  }
  memset(st, 0, num_queries * sizeof(int));
  return attack0(r_filt, t_filts, n_filts, st, filter_bytes, num_queries,
                 n, 1, 0, dbg, test_ct);
}

int test(int k, int m, int n, int r, int q, int dbg, int seed) {
  srand(seed);
  int filter_bytes = (1<<m)/8;
  int num_queries = 1<<q;

  char *r_filt = malloc(filter_bytes * sizeof(char));

  char **t_filts = malloc(num_queries * sizeof(char*));
  for (int i = 0; i < num_queries; i++) {
    t_filts[i] = malloc(filter_bytes * sizeof(char));
  }

  char **n_filts = malloc(num_queries * sizeof(char*));
  for (int i = 0; i < num_queries; i++) {
    n_filts[i] = malloc(filter_bytes * sizeof(char));
  }
  int *st = malloc(num_queries * sizeof(int));

  // Sample the target filter.
  random_filter(r_filt, k, m, r);
  if (dbg) {
    printf("R= "); print_filter(r_filt, m);
  }

  // Sample the test filters.
  for (int i = 0; i < num_queries; i++) {
    random_filter(t_filts[i], k, m, 1);
    if (dbg) {
      printf("M%d=", i); print_filter(t_filts[i], m);
    }
  }

  //Run attacks.
  int res = attack(r_filt, t_filts, n_filts, st, k, m, n, r, q, dbg, NULL);
  int ref_res = ref_attack(r_filt, t_filts, k, m, n, r, q, 1000);

  free(r_filt);
  for (int i = 0; i < num_queries; i++) {
    free(t_filts[i]);
    free(n_filts[i]);
  }
  free(t_filts);
  free(n_filts);
  free(st);

  if (!res && ref_res) {
    fprintf(stderr, "test fails: seed=%d\n", seed);
    return 0;
  }
  return 1;
}

int multi_test(int k, int m, int n, int r, int q, int trials) {
  int bad = 0;
  for (int i = 0; i < trials; i++) {
    if (!test(k, m, n, r, q, 0, rand())) {
      bad++;
    }
  }
  if (bad == 0) {
    printf("pass\n");
    return 1;
  } else {
    fprintf(stderr, "%d of %d failed\n", bad, trials);
    return 0;
  }
}

double sim(int k, int m, int n, int r, int q, int trials, double *avg_tests) {
  int filter_bytes = (1<<m)/8;
  int num_queries = 1<<q;

  char *r_filt = malloc(filter_bytes * sizeof(char));

  char **t_filts = malloc(num_queries * sizeof(char*));
  for (int i = 0; i < num_queries; i++) {
    t_filts[i] = malloc(filter_bytes * sizeof(char));
  }

  char **n_filts = malloc(num_queries * sizeof(char*));
  for (int i = 0; i < num_queries; i++) {
    n_filts[i] = malloc(filter_bytes * sizeof(char));
  }
  int *st = malloc(num_queries * sizeof(int));

  int good = 0, ct = 0;
  for (int i = 0; i < trials; i++) {
    // Sample the target filter.
    random_filter(r_filt, k, m, r);

    // Sample the test filters.
    for (int i = 0; i < num_queries; i++) {
      random_filter(t_filts[i], k, m, 1);
    }

    // Execute the attack.
    good += attack(r_filt, t_filts, n_filts, st, k, m, n, r, q, 0, &ct);
  }

  *avg_tests = (double) ct / trials;

  free(r_filt);
  for (int i = 0; i < num_queries; i++) {
    free(t_filts[i]);
    free(n_filts[i]);
  }
  free(t_filts);
  free(n_filts);
  free(st);

  return (double) good / trials;
}

int main(int argc, const char **argv) {

  // Seed the RNG. (This suffices for simulating the attack.)
  int seed = time(NULL);
  srand(seed);

  // Default parameters.
  int trials = 1000,
      k=4,   // number of hashes
      m=10,  // 2^m = filter bits
      n=100, // set size
      r=1,   // error parameter
      q=10;  // 2^q = number of queries

  if ((1<<q) < n) {
    fprintf(stderr, "Bad parameters.\n");
    return 0;
  }

  // Run tests.
  if (!multi_test(2,4,2,1,5,1000)) {
    return 1;
  }

  // Run simulations.
  double adv, avg_tests;

  printf("\nvarying, adv, min_adv, avg_tests\n");
  for (int vk=1; vk <= 8; vk += 1) {
    adv = sim(vk, m, n, r, q, trials, &avg_tests);
    printf(" k=%-5d %-.3f %-.2f\n",
        vk, adv, avg_tests);
  }

  printf("\n");
  for (int vm = 9; vm <= 15; vm += 1) {
    adv = sim(k, vm, n, r, q, trials, &avg_tests);
    printf(" m=%-5d %-.3f %-.2f\n",
        1<<vm, adv, avg_tests);
  }

  printf("\n");
  for (int vn = k*r; vn <= (k*r)+5; vn += 1) {
    adv = sim(k, m, vn, r, q, trials, &avg_tests);
    printf(" n=%-5d %-.3f %-.2f\n",
        vn, adv, avg_tests);
  }

  printf("\n");
  for (int vr = r; vr <= 7; vr += 1) {
    adv = sim(k, m, n, vr, q, trials, &avg_tests);
    printf(" r=%-5d %-.3f %-.2f \n",
        vr, adv, avg_tests);
  }

  printf("\n");
  for (int vq = 7; vq <= 13; vq += 1) {
    adv = sim(k, m, n, r, vq, trials, &avg_tests);
    printf(" q=%-5d %-.3f %-.2f \n",
        1<<vq, adv, avg_tests);
  }

  return 0;
}
