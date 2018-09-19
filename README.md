# Introduction

This repository contains our computational experiment for testing various instances of density matrix and measurements for violation of C2(b) coherence property. We tested different kinds of norms and this program was used as tool for checking for violations of C2(b).

## Get started

Install all the requirements using

`pip install -r requirements.txt`

From there you can perform three kinds of computational experiment. Random sampling, test using precomputed data and meta-heuristic test. Each of those has a dedicated section below.

## Test using random samples

Numerical test of `C2(b)` with `schatten-1` norm and `L1` norm (for debugging) on one million instances.

Run

`python rsample.py`

Expected result

```
L = 2, N = 2
100%|███████████████| 100000/100000 [07:10<00:00, 232.41it/s]

total violations so far: 0

L = 2, N = 4
100%|███████████████| 100000/100000 [08:38<00:00, 192.94it/s]

total violations so far: 0

L = 4, N = 2
100%|██████████████████████████████████| 100000/100000 [14:13<00:00, 117.11it/s]

total violations so far: 0

L = 4, N = 4
100%|███████████████████████████████████| 100000/100000 [18:10<00:00, 91.72it/s]

total violations so far: 0

L = 8, N = 2
100%|███████████████████████████████████| 100000/100000 [30:34<00:00, 54.51it/s]

total violations so far: 0

L = 8, N = 4
100%|███████████████████████████████████| 100000/100000 [42:40<00:00, 39.06it/s]

total violations so far: 0

L = 16, N = 2
100%|█████████████████████████████████| 100000/100000 [1:19:11<00:00, 21.05it/s]

total violations so far: 0

L = 16, N = 4
100%|█████████████████████████████████| 100000/100000 [1:59:23<00:00, 13.96it/s]

total violations so far: 0

L = 32, N = 2
100%|█████████████████████████████████| 100000/100000 [3:53:51<00:00,  7.13it/s]

total violations so far: 0

L = 32, N = 4
100%|█████████████████████████████████| 100000/100000 [6:16:29<00:00,  4.43it/s]

total violations so far: 0


total number of instances: 1000000


total number of violated instances: 0


number of violating instances per test:

+--- norm ----+-- Pn / P --+-- B + iBd --+
| L01         |          0 |           0 |
| Schatten-01 |          0 |           0 |
+-------------+------------+-------------+

```

## Test using pre-generated problem instances

At the moment of development of this computational experiment [qutip](http://qutip.org) library doesn't give control of randomness, we cannot pass our own seed as parameter thus it is harder to reproduce results. As a dirty solution we implemented a script that pre-generates problem instances and stores them in files. You can feel free to generate as many as you need, during the test instances are streamed from files avoiding wasting memory.

We pre-generated some instances, look up [instances directory](https://github.com/nyuqtl/schatten/tree/master/instances).

### Generating instances

Generate set of random density matrices

`python density_gen.py`

Generate random measurements

`python x_gen.py`

### Running test on pre-generated instances

Once set of instances is generated you can run the test script

`python fsample.py`

Output is compatible with one from random sampling

## Simulated annealing

In order to avoid relying entirely on random sampling we also implemented a meta-heuristic approach based on simulated annealing.

The way it works is initially it random samples for density matrix for which C2(b) is at equality, and then it uses simulated annealing to explore local neighborhood of solutions to search for violation. Simulated annealing part is based on [simanneal module](https://github.com/perrygeo/simanneal/)

Reason why we decided to implement it that way was simulated annealing performed well locally, solution was improving up to 30%, however local search was not aggressive enough to lift C2(b) to equality, even with extreme initial temperatures.

Run

`python anneal.py`

Output

```
random sampling for density matrix

0.0444724779802278
0.036675873109657975
0.03453342161596007
0.03423210669316901
0.025258985656120198
0.024990294827152662
0.009029733290144004
0.0070368609534336435
0.004431255920851647
0.003605852687035574
0.0031558219717035575
0.0

found tightest C2(b) bound which is approximately 0.
beginning simulated annealing to find measurement violating C2(b) as local optimum

 Temperature        Energy    Accept   Improve     Elapsed   Remaining
     2.50000          0.00   100.00%     0.00%     0:02:07     0:00:00


no C2(b) violation found
```
