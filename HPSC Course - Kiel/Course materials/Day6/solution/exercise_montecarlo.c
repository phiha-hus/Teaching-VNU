

#include "basic.h"

static const int rng_m = (1u << 31) - 1;
static const int rng_a = 16807;
static const int rng_seed = 58854338;
static const int rng_b = 0;
static const float invRAND_MAX = 1.0f / ((1u << 31) - 1);

/* very simple pseudo-random-number-generator */
float myrng(int *s)
{
  *s = (rng_a * *s + rng_b) % rng_m;
  return (float)*s * invRAND_MAX;
}

/* monte carlo method for computing the value of pi, single threaded */
int montecarlo_pi(int samples)
{
  float x, y;
  int i;
  int seed = rng_seed;

  int hit;

  hit = 0;
  for (i = 0; i < samples; ++i)
  {
    x = myrng(&seed);
    y = myrng(&seed);

    if (x * x + y * y <= 1)
    {
      hit++;
    }
  }

  return hit;
}

int offset = (64 / sizeof(int));

/* monte carlo method for computing the value of pi, multi threaded,
 the good one */
int montecarlo_pi_omp(int samples)
{
  int i;
  int seeds[8 * offset];
  int *p;

  int hit;

  hit = 0;
#pragma omp parallel
  {
    float x, y;
    seeds[omp_get_thread_num() * offset] =
        (omp_get_thread_num() + 1) * rng_seed;

#pragma omp for reduction(+ \
                          : hit)
    for (i = 0; i < samples; ++i)
    {
      p = seeds + omp_get_thread_num() * offset;
      x = myrng(p);
      y = myrng(p);

      if (x * x + y * y <= 1.0f)
      {
        hit++;
      }
    }
  }

  return hit;
}

/* monte carlo method for computing the value of pi, multi threaded,
 the bad one */
int montecarlo_pi_omp_FS(int samples)
{
  float x[2];
  int i;
  int seeds[8];

  int hit;

  hit = 0;
#pragma omp parallel num_threads(8)
  {
    // possible False Sharing
    seeds[omp_get_thread_num()] = (omp_get_thread_num() + 1) * rng_seed;

#pragma omp for
    for (i = 0; i < samples; ++i)
    {
      // possible False Sharing
      x[0] = myrng(seeds + omp_get_thread_num());
      x[1] = myrng(seeds + omp_get_thread_num());

      if (x[0] * x[0] + x[1] * x[1] <= 1.0f)
      {
        // Atomic reduces parallel efficiency,
        // better use private variables and combine results at the end
#pragma omp atomic
        hit++;
      }
    }
  }

  return hit;
}

int main(int argc, char **argv)
{
  int samples;
  float apprx_pi, apprx_pi_omp, exact_pi;
  pstopwatch sw;
  double t;

  if (argc != 2)
  {
    fprintf(stderr, " usage: montecarlo <samples>\n");
  }
  samples = atoi(argv[1]);

  sw = new_stopwatch();

  exact_pi = M_PI;

  printf("Computing pi sequential\n");
  start_stopwatch(sw);
  apprx_pi = 4.0f * (float)montecarlo_pi(samples) / (float)samples;
  t = stop_stopwatch(sw);
  printf("  %.2f s\n", t);

  printf("Computing pi parallel\n");
  start_stopwatch(sw);
  apprx_pi_omp = 4.0f * (float)montecarlo_pi_omp_FS(samples) / (float)samples;
  t = stop_stopwatch(sw);
  printf("  %.2f s\n", t);

  printf("exact value of pi:  %.15f\n", exact_pi);
  printf("apprx1 value of pi: %.15f\n", apprx_pi);
  printf("relative error:    %.5e\n", fabsf(exact_pi - apprx_pi) / exact_pi);
  printf("apprx2 value of pi: %.15f\n", apprx_pi_omp);
  printf("relative error:    %.5e\n",
         fabsf(exact_pi - apprx_pi_omp) / exact_pi);

  printf("Computing pi parallel\n");
  start_stopwatch(sw);
  apprx_pi_omp = 4.0f * (float)montecarlo_pi_omp(samples) / (float)samples;
  t = stop_stopwatch(sw);
  printf("  %.2f s\n", t);

  printf("apprx3 value of pi: %.15f\n", apprx_pi_omp);
  printf("relative error:    %.5e\n",
         fabsf(exact_pi - apprx_pi_omp) / exact_pi);

  return 0;
}
