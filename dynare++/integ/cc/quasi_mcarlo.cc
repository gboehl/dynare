// Copyright 2005, Ondra Kamenik

#include "quasi_mcarlo.hh"

#include <cmath>

/* Here in the constructor, we have to calculate a maximum length of
   |coeff| array for a given |base| and given maximum |maxn|. After
   allocation, we calculate the coefficients. */

RadicalInverse::RadicalInverse(int n, int b, int mxn)
  : num(n), base(b), maxn(mxn),
    coeff((int) (floor(log((double)maxn)/log((double)b))+2), 0)
{
  int nr = num;
  j = -1;
  do
    {
      j++;
      coeff[j] = nr % base;
      nr = nr / base;
    }
  while (nr > 0);
}

/* This evaluates the radical inverse. If there was no permutation, we have to calculate
   $$
   {c_0\over b}+{c_1\over b^2}+\ldots+{c_j\over b^{j+1}}
   $$
   which is evaluated as
   $$
   \left(\ldots\left(\left({c_j\over b}\cdot{1\over b}+{c_{j-1}\over b}\right)
   \cdot{1\over b}+{c_{j-2}\over b}\right)
   \ldots\right)\cdot{1\over b}+{c_0\over b}
   $$
   Now with permutation $\pi$, we have
   $$
   \left(\ldots\left(\left({\pi(c_j)\over b}\cdot{1\over b}+
   {\pi(c_{j-1})\over b}\right)\cdot{1\over b}+
   {\pi(c_{j-2})\over b}\right)
   \ldots\right)\cdot{1\over b}+{\pi(c_0)\over b}
   $$ */

double
RadicalInverse::eval(const PermutationScheme &p) const
{
  double res = 0;
  for (int i = j; i >= 0; i--)
    {
      int cper = p.permute(i, base, coeff[i]);
      res = (cper + res)/base;
    }
  return res;
}

/* We just add 1 to the lowest coefficient and check for overflow with respect
   to the base. */
void
RadicalInverse::increase()
{
  // todo: raise if |num+1 > maxn|
  num++;
  int i = 0;
  coeff[i]++;
  while (coeff[i] == base)
    {
      coeff[i] = 0;
      coeff[++i]++;
    }
  if (i > j)
    j = i;
}

/* Debug print. */
void
RadicalInverse::print() const
{
  printf("n=%d b=%d c=", num, base);
  coeff.print();
}

/* Here we have the first 170 primes. This means that we are not able
   to integrate dimensions greater than 170. */

int HaltonSequence::num_primes = 170;
int HaltonSequence::primes[] = {
  2,     3,     5,     7,    11,    13,    17,    19,    23,    29,
  31,    37,    41,    43,    47,    53,    59,    61,    67,    71,
  73,    79,    83,    89,    97,   101,   103,   107,   109,   113,
  127,   131,   137,   139,   149,   151,   157,   163,   167,   173,
  179,   181,   191,   193,   197,   199,   211,   223,   227,   229,
  233,   239,   241,   251,   257,   263,   269,   271,   277,   281,
  283,   293,   307,   311,   313,   317,   331,   337,   347,   349,
  353,   359,   367,   373,   379,   383,   389,   397,   401,   409,
  419,   421,   431,   433,   439,   443,   449,   457,   461,   463,
  467,   479,   487,   491,   499,   503,   509,   521,   523,   541,
  547,   557,   563,   569,   571,   577,   587,   593,   599,   601,
  607,   613,   617,   619,   631,   641,   643,   647,   653,   659,
  661,   673,   677,   683,   691,   701,   709,   719,   727,   733,
  739,   743,   751,   757,   761,   769,   773,   787,   797,   809,
  811,   821,   823,   827,   829,   839,   853,   857,   859,   863,
  877,   881,   883,   887,   907,   911,   919,   929,   937,   941,
  947,   953,   967,   971,   977,   983,   991,   997,  1009,  1013
};

/* This takes first |dim| primes and constructs |dim| radical inverses
   and calls |eval|. */

HaltonSequence::HaltonSequence(int n, int mxn, int dim, const PermutationScheme &p)
  : num(n), maxn(mxn), per(p), pt(dim)
{
  // todo: raise if |dim > num_primes|
  // todo: raise if |n > mxn|
  for (int i = 0; i < dim; i++)
    ri.emplace_back(num, primes[i], maxn);
  eval();
}

const HaltonSequence &
HaltonSequence::operator=(const HaltonSequence &hs)
{
  num = hs.num;
  maxn = hs.maxn;
  ri.clear();
  for (const auto & i : hs.ri)
    ri.emplace_back(i);
  pt = hs.pt;
  return *this;
}

/* This calls |RadicalInverse::increase| for all radical inverses and
   calls |eval|. */

void
HaltonSequence::increase()
{
  for (auto & i : ri)
    i.increase();
  num++;
  if (num <= maxn)
    eval();
}

/* This sets point |pt| to radical inverse evaluations in each dimension. */
void
HaltonSequence::eval()
{
  for (unsigned int i = 0; i < ri.size(); i++)
    pt[i] = ri[i].eval(per);
}

/* Debug print. */
void
HaltonSequence::print() const
{
  for (const auto & i : ri)
    i.print();
  printf("point=[ ");
  for (unsigned int i = 0; i < ri.size(); i++)
    printf("%7.6f ", pt[i]);
  printf("]\n");
}

qmcpit::qmcpit()
  : spec(NULL), halton(NULL), sig(NULL)
{
}

qmcpit::qmcpit(const QMCSpecification &s, int n)
  : spec(&s), halton(new HaltonSequence(n, s.level(), s.dimen(), s.getPerScheme())),
    sig(new ParameterSignal(s.dimen()))
{
}

qmcpit::qmcpit(const qmcpit &qpit)
  : spec(qpit.spec), halton(NULL), sig(NULL)
{
  if (qpit.halton)
    halton = new HaltonSequence(*(qpit.halton));
  if (qpit.sig)
    sig = new ParameterSignal(qpit.spec->dimen());
}

qmcpit::~qmcpit()
{
  if (halton)
    delete halton;
  if (sig)
    delete sig;
}

bool
qmcpit::operator==(const qmcpit &qpit) const
{
  return (spec == qpit.spec)
    && ((halton == NULL && qpit.halton == NULL)
        || (halton != NULL && qpit.halton != NULL && halton->getNum() == qpit.halton->getNum()));
}

const qmcpit &
qmcpit::operator=(const qmcpit &qpit)
{
  spec = qpit.spec;
  if (halton)
    delete halton;
  if (qpit.halton)
    halton = new HaltonSequence(*(qpit.halton));
  else
    halton = NULL;
  return *this;
}

qmcpit &
qmcpit::operator++()
{
  // todo: raise if |halton == null || qmcq == NULL|
  halton->increase();
  return *this;
}

double
qmcpit::weight() const
{
  return 1.0/spec->level();
}

qmcnpit::qmcnpit()
  : qmcpit(), pnt(NULL)
{
}

qmcnpit::qmcnpit(const QMCSpecification &s, int n)
  : qmcpit(s, n), pnt(new Vector(s.dimen()))
{
}

qmcnpit::qmcnpit(const qmcnpit &qpit)
  : qmcpit(qpit), pnt(NULL)
{
  if (qpit.pnt)
    pnt = new Vector(*(qpit.pnt));
}

qmcnpit::~qmcnpit()
{
  if (pnt)
    delete pnt;
}

const qmcnpit &
qmcnpit::operator=(const qmcnpit &qpit)
{
  qmcpit::operator=(qpit);
  if (pnt)
    delete pnt;
  if (qpit.pnt)
    pnt = new Vector(*(qpit.pnt));
  else
    pnt = NULL;
  return *this;
}

/* Here we inccrease a point in Halton sequence ant then store images
   of the points in |NormalICDF| function. */

qmcnpit &
qmcnpit::operator++()
{
  qmcpit::operator++();
  for (int i = 0; i < halton->point().length(); i++)
    (*pnt)[i] = NormalICDF::get(halton->point()[i]);
  return *this;
}

/* Clear from code. */
int
WarnockPerScheme::permute(int i, int base, int c) const
{
  return (c+i) % base;
}

/* Clear from code. */
int
ReversePerScheme::permute(int i, int base, int c) const
{
  return (base-c) % base;
}
