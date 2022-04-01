CX4220/CSE 6220 Programming Assignment 2
========================================

**DO NOT FORK THIS REPO. DO NOT PUSH PROJECT SOLUTIONS PUBLICLY.
THIS IS AN ACADEMIC INTEGRITY VIOLATION.**

## Code structure

All the code is located at the root level of the project.

There are multiple header and .cpp files, your implementation will go
into `collective.cpp`. You will only submit `collective.cpp`.

- `collective.h`: This header file contains the function declarations for `HPC_Bcast()` and `HPC_Prefix()`.
- `Makefile`: contains commands for compiling and running the tests.
- `mpi_tests.cpp`: This holds tests that will be used for the autograder. You
  can look over this file to see how your functions will be called. You can
  modify this file to run your own tests or to help debug your code however you
  want, but you will not submit this file.
- `utils.h`,`utils.cpp`: These contain a function for printing out integer
  arrays in a somewhat ordered fashion. It may be helpful for debugging.
- `example_pbs_script.pbs`: This is a simple PBS script for running the tests.
  You can use it as a starting point for running your own tests on PACE.
- `gtest/`,`gtest-mpi-listener/`,`mpi_gtest.cpp`: These are just libraries for
  the test framework. You can ignore them completely.

## Compiling and Testing

In order to compile the tests in `mpi_tests`, simply run `make`.

You can run the tests using `make test_all`, `make test_broadcast`, or `make test_prefix`.

For example,
```
    >$ make test_broadcast 
mpic++ -Wall --std=c++17 -g -O3 -I"./gtest" -I"./gtest/src" -I"./gtest/include" -c mpi_gtest.cpp
mpic++ -Wall --std=c++17 -g -O3 -I"./gtest" -I"./gtest/src" -I"./gtest/include"  -o mpi_tests mpi_tests.o mpi_gtest.o gtest-all.o collective.o utils.o
mpirun -np 4 ./mpi_tests --gtest_filter='*Bcast*'
Note: Google Test filter = *Bcast*
[==========] Running 8 tests from 2 test suites.
[----------] Global test environment set-up.
[----------] 4 tests from BcastInt_ZeroRoot/BcastTest
[ RUN      ] BcastInt_ZeroRoot/BcastTest.BcastN/0
[       OK ] BcastInt_ZeroRoot/BcastTest.BcastN/0 (0 ms)
[ RUN      ] BcastInt_ZeroRoot/BcastTest.BcastN/1
[       OK ] BcastInt_ZeroRoot/BcastTest.BcastN/1 (0 ms)
[ RUN      ] BcastInt_ZeroRoot/BcastTest.BcastN/2
[       OK ] BcastInt_ZeroRoot/BcastTest.BcastN/2 (0 ms)
[ RUN      ] BcastInt_ZeroRoot/BcastTest.BcastN/3
[       OK ] BcastInt_ZeroRoot/BcastTest.BcastN/3 (1 ms)
[----------] 4 tests from BcastInt_ZeroRoot/BcastTest (1 ms total)

[----------] 4 tests from BcastInt_NonzeroRoot/BcastTest
[ RUN      ] BcastInt_NonzeroRoot/BcastTest.BcastN/0
[       OK ] BcastInt_NonzeroRoot/BcastTest.BcastN/0 (0 ms)
[ RUN      ] BcastInt_NonzeroRoot/BcastTest.BcastN/1
[       OK ] BcastInt_NonzeroRoot/BcastTest.BcastN/1 (0 ms)
[ RUN      ] BcastInt_NonzeroRoot/BcastTest.BcastN/2
[       OK ] BcastInt_NonzeroRoot/BcastTest.BcastN/2 (0 ms)
[ RUN      ] BcastInt_NonzeroRoot/BcastTest.BcastN/3
[  SKIPPED ] BcastInt_NonzeroRoot/BcastTest.BcastN/3 (0 ms)
[----------] 4 tests from BcastInt_NonzeroRoot/BcastTest (0 ms total)

[----------] Global test environment tear-down
Finalizing MPI...
[==========] 8 tests from 2 test suites ran. (50 ms total)
[  PASSED  ] 7 tests.
[  SKIPPED ] 1 test, listed below:
[  SKIPPED ] BcastInt_NonzeroRoot/BcastTest.BcastN/3
```

To have more control over which tests are run, you can run `./mpi_tests` with
gtest arguments. Run `./mpi_tests --help` to see the list of possible
commandline arguments. For example, `--gtest_filter` can be used to run a
specific subset of the tests.


```
   >$ mpirun -np 4 ./mpi_tests --gtest_filter='PrefixSumTestSuite/PrefixSumTest.Add/2'
Note: Google Test filter = PrefixSumTestSuite/PrefixSumTest.Add/2
[==========] Running 1 test from 1 test suite.
[----------] Global test environment set-up.
[----------] 1 test from PrefixSumTestSuite/PrefixSumTest
[ RUN      ] PrefixSumTestSuite/PrefixSumTest.Add/2
mpi_tests.cpp:114: Failure
Failed
[Rank 1/4] Expected equality of these values:
[Rank 1/4]   recvbuffer[i-rank*n]
[Rank 1/4]     Which is: 43
[Rank 1/4]   i*(i+1)/2
[Rank 1/4]     Which is: 946
[Rank 1/4]  recvbuffer[0] is wrong on rank 1

mpi_tests.cpp:114: Failure
Failed
[Rank 2/4] Expected equality of these values:
[Rank 2/4]   recvbuffer[i-rank*n]
[Rank 2/4]     Which is: 86
[Rank 2/4]   i*(i+1)/2
[Rank 2/4]     Which is: 3741
[Rank 2/4]  recvbuffer[0] is wrong on rank 2

mpi_tests.cpp:114: Failure
Failed
[Rank 3/4] Expected equality of these values:
[Rank 3/4]   recvbuffer[i-rank*n]
[Rank 3/4]     Which is: 129
[Rank 3/4]   i*(i+1)/2
[Rank 3/4]     Which is: 8385
[Rank 3/4]  recvbuffer[0] is wrong on rank 3

[  FAILED  ] PrefixSumTestSuite/PrefixSumTest.Add/2, where GetParam() = (3, 43) (0 ms)
[----------] 1 test from PrefixSumTestSuite/PrefixSumTest (0 ms total)

[----------] Global test environment tear-down
Finalizing MPI...
[==========] 1 test from 1 test suite ran. (51 ms total)
[  PASSED  ] 0 tests.
[  FAILED  ] 1 test, listed below:
[  FAILED  ] PrefixSumTestSuite/PrefixSumTest.Add/2, where GetParam() = (3, 43)

 1 FAILED TEST
```

