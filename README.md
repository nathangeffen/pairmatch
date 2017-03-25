# Pair matching algorithms

This is a PhD research project.

The code is in C++ and has been tested on GNU/Linux machines using g++ (version
5 and up should be fine).

To approximately replicate the data in submitted paper, execute
run_tests.sh. Exact replication isn't likely because of stochastic effects.

To execute the main program, go to the code directory and run *make release* or
*make release_attract* (for the two different agent based models described in
the paper).

## Command line options

- -n <unsigned int> : Number of agents
- -k <unsigned int> : Value of k in several algorithms
- -c <unsigned int> : Number of clusters in CSPM
- -y <unsigned int> : Number of age years to look forward and back in DCPM
- -i <unsigned int> : Number of iterations per simulation
- -r <unsigned int> : Number of simulations per algorithm
- -s <unsigned int> : Random number seed
- -f <unsigned int> : Identifying number in output of the first simulation
- -a [R|N|W|D|C|B]  : Algorithms to run: R=RPM, N=RKPM, W=WSPM, D=DCPM, C=CSPM,
- -B=BFPM
- -A <float [0..1]> : Attractor parameter for ATTRACT_REJECT
- -R <float [0..1]> : Rejector parameter for ATTRACT_REJECT
- -o <STRING>       : Identifying string for output (see run_tests.sh)
- -d <STRING>       : Identifying string for output (see run_tests.sh)

Also see:
http://nathangeffen.webfactional.com/partnermatching/partnermatching.html
