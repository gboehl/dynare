// Copyright 2007, Ondra Kamenik

#include "random.hh"

#include <cstdlib>
#include <limits>
#include <cmath>

SystemRandomGenerator system_random_generator;

int
RandomGenerator::int_uniform()
{
  double s = std::numeric_limits<int>::max()*uniform();
  return (int) s;
}

/* This implements Marsaglia Polar Method. */
double
RandomGenerator::normal()
{
  double x1, x2;
  double w;
  do
    {
      x1 = 2*uniform()-1;
      x2 = 2*uniform()-1;
      w = x1*x1 + x2*x2;
    }
  while (w >= 1.0 || w < 1.0e-30);
  return x1*std::sqrt((-2.0*std::log(w))/w);
}

double
SystemRandomGenerator::uniform()
{
#if !defined(__MINGW32__)
  return drand48();
#else
  return ((double) rand())/RAND_MAX;
#endif
}

void
SystemRandomGenerator::initSeed(int seed)
{
#if !defined(__MINGW32__)
  srand48(seed);
#else
  srand(seed);
#endif
}
