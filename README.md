# heptools
High energy physics utilities
This package provides
* unit conversion
* beta functions for SM
* group manipulations
and so on. (TBD)

# samples
## unit calculation

```python:
from heptools.units import *

# 1.52e+15 * GeV**(-1)
(30 * cm).in_GeV()

# 6.58e-16 * GeV
(30 * cm).in_GeV().inverse()

# 6.58e-25 * s
(GeV**-1).in_s()

# 1.00e-39 * cm**2
fb.in_cm()

# 3.89e+11 * GeV**2
fb.inverse().in_GeV()
```

## Lie algebra
We should improve the algorithm referring to 1206.6379.
Also, we are going to implement representation classes.

```python:
from heptools.lie import *

su3 = su(3)
su3.index([3,0])

```

