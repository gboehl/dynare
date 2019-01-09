// Copyright 2007, Ondra Kamenik

// Random number generation

#ifndef RANDOM_H
#define RANDOM_H

/* This is a general interface to an object able to generate random
   numbers. Subclass needs to implement |uniform| method, other is, by
   default, implemented here. */
class RandomGenerator
{
public:
  virtual double uniform() = 0;
  int int_uniform();
  double normal();
};

/* This implements |RandomGenerator| interface with system |drand| or
   |rand|. It is not thread aware. */
class SystemRandomGenerator : public RandomGenerator
{
public:
  double uniform() override;
  void initSeed(int seed);
};

extern SystemRandomGenerator system_random_generator;

#endif
