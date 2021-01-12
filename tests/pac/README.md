Collection of `*.mod` files for testing the PAC routines in Dynare. To run all the tests sequentially and check that all the `*.mod` files pass the tests, just use the matlab function in the base directory:

```matlab
>> run_all_tests()
```

If all goes well, the output should terminate with something like:

```example
Testsuite results (PAC model):

var-1		             PASS (2.0479s)
var-2		             PASS (1.9601s)
var-3		             PASS (1.9826s)
var-4		             PASS (2.0079s)
trend-component-1		 PASS (2.2214s)
trend-component-2		 PASS (2.2195s)
trend-component-3		 PASS (2.3003s)
trend-component-4		 PASS (10.5143s)
trend-component-5		 PASS (2.1538s)
trend-component-6		 PASS (2.4203s)
trend-component-7		 PASS (2.7112s)
trend-component-9		 PASS (2.2164s)
trend-component-10		 PASS (2.2593s)
trend-component-11		 PASS (0.69409s)
trend-component-12		 PASS (2.2663s)
trend-component-13a		 PASS (0.4119s)
trend-component-13b		 PASS (0.39554s)
```
