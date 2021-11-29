# snowfake

[![PyPI status](https://img.shields.io/pypi/status/snowfake.svg)](https://pypi.org/project/snowfake//)
[![PyPI versions](https://img.shields.io/pypi/pyversions/snowfake.svg)](https://pypi.org/project/snowfake//)
[![PyPI license](https://img.shields.io/pypi/l/snowfake.svg)](https://pypi.org/project/snowfake/)

Make Gravner-Griffeath "snowfakes"! This code implements:

> Janko Gravner, David Griffeath (2008). Modeling snow crystal growth II: A mesoscopic lattice map with plausible dynamics. _Physica D: Nonlinear Phenomena_ **237** (3), p 385-404. [DOI: 10.1016/j.physd.2007.09.008](https://doi.org/10.1016/j.physd.2007.09.008).


## Installation

You can install this package with `pip` (be careful not to type "snowflake"):

    pip install snowfake

#### Experimental

Installing `scikit-image` allows you to use a different affine transformation, but I haven't figured out yet if it's better or not. 

    pip install snowfake[skimage]


## Example

You can produce a randome snowfake with:

    import snowfake
    s = snowfake.random()

This code produces the crystal in Figure 5b of the Gravner & Griffeath (2008):

```python
from snowfake import Snowfake

params =  {
    'ρ': 0.35,  # or 'rho': 0.35 if you prefer...
    'β': 1.4,
    'α': 0.001,
    'θ': 0.015,
    'κ': 0.05,
    'μ': 0.015,
    'γ': 0.01,
    'σ': 0.00005,
    'random': False,
}
s = Snowfake(size=801, **params)
s.grow()
s.plot()
```

The various physical parameter arrays are available as `s.a` (attachment flag), `s.b` (boundary mass), `s.c` (the crystal itself) and `s.d` (the vapour). The arrays exist on hexgrids; you can rectify them with, for example, `s.rectify('c')`.

The parameter `σ` (note that you can also spell out `sigma` if you prefer) can be a 1D array with one sample per epoch. This will vary the vapour density `ρ` through _time_. The parameter `ρ` can be a 2D array of shape `(size, size)`; this will vary the initial vapour density through _space_.


## Testing

You can run the tests (requires `pytest` and `pytest-cov`) with 

    python run_tests.py


---

&copy; 2021 Agile Scientific, openly licenced under Apache 2.0