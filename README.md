# heputils
High energy physics utilities
This package provides
* unit conversion
* beta functions for SM
* group manipulations
and so on. (TBD)

# samples
## unit calculation

```python:
from heputils.units import *

# 1.52e+15 * GeV**(-1)
(30 * cm).in_(GeV)

# 6.58e-16 * GeV
(30 * cm).in_(GeV).inverse()

# 6.58e-25 * s
(GeV**-1).in_(s)
```

