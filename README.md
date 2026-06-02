# Factor-Revealing LPs for the Unsplittable CVRP

## Overview

* `simple_lp.py`: Implements the factor-revealing LP for the simplified framework (Theorem 3).

* `final_lp.py`: Implements the factor-revealing LPs for the full framework (Theorems 4 and 5).

* `results/`: Gurobi logs for the LPs in the full framework, corresponding to the results reported in Tables 3--6.

---

## Python Files

### `simple_lp.py`

**Input:** `<N>`

* `N`: Demand discretization parameter $N$ used in the paper. It should be a multiple of $12$.

**Output:** `<N, value, runtime>`

* LP value and running time.

**Note:** $M=1000$ by default.

### `final_lp.py`

**Input:** `<metricflag, randomflag, generalflag, N>`

* `metricflag`: `1` = Metric case, `0` = Asymmetric case.
* `randomflag`: `1` = Randomized algorithm, `0` = Deterministic algorithm.
* `generalflag`: `1` = General-capacity case, `0` = Fixed-capacity case.
* `N`: Demand discretization parameter $N$ used in the paper. It should be a multiple of $12$.

**Output:** `<N, value, runtime>`

* LP value and running time.

**Note:** $N_{\max}=4$ and $M=1000$ by default.

For example:

* `<metricflag, randomflag, generalflag> = <0,1,1>` corresponds to the LP of `R-General-Alg` in the asymmetric case.
* `<metricflag, randomflag, generalflag> = <0,1,0>` corresponds to the LP of `R-Fixed-Alg` in the asymmetric case.
* Other combinations are interpreted analogously.

---

## Computational Environment

The LPs were implemented in Python 3.8.10 and solved using Gurobi Optimizer 9.5.1.

All experiments were conducted on a server with:

* Intel Xeon Gold 6226R CPU @ 2.90 GHz
* 512 GB RAM
* Ubuntu Server 20.04 LTS

---

## Gurobi Parameters

To improve numerical stability when solving large-scale LPs, we use Gurobi's barrier (interior-point) method together with several numerical-stability settings:

```python
m.setParam("Threads", 4)
m.setParam("Method", 2)
m.setParam("BarOrder", 1)
m.setParam("BarHomogeneous", 1)
m.setParam("NumericFocus", 2)
m.setParam("Crossover", 0)
```

---

## Requirements

* Python 3.8+ (tested with Python 3.8.10)
* Gurobi Optimizer (tested with version 9.5.1)
* `gurobipy` (with a valid Gurobi license)

---

## How to Run

```bash
python simple_lp.py <N>

python final_lp.py <metricflag> <randomflag> <generalflag> <N>
```

### Examples

Input:

```bash
python simple_lp.py 60
```

Output:

```text
N: 60, value: 1.6653220916579425, runtime: 0.3686387538909912 seconds
```

Input:

```bash
python final_lp.py 0 1 1 60
```

Output:

```text
N: 60, value: 1.6501091849268639, runtime: 44.63564705848694 seconds
```
