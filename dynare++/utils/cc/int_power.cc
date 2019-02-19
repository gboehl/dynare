int
power(int a, int b)
{
  int res = 1;
  while (b > 0)
    {
      if (b & 1)
        res *= a;
      a *= a;
      b >>= 1;
    }
  return res;
}
