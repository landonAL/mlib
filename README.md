# MLib

> [!WARNING]
> This library is a work in progress, and anything is subject to change without notice. Use this library at your own risk!

A mathematical library that does not have any dependencies and aims to provide advanced mathematical functionality.

The library itself does not need to be built. You can simply copy-paste [./mlib.h](./mlib.h) to your project and `#include` it. (Because the truly reusable code is the one that you can simply copy-paste).

MLib is a classical [stb-style](https://github.com/nothings/stb) single header library. That is by default it acts like a header, but if you `#define MLIB_IMPLEMENTATION` prior to including the library it acts like a source file.

## Quick Example

```c
#include <stdio.h>

#define MLIB_IMPLEMENTATION
#include "mlib.h"

int main(void)
{
    for (int i = 1; i < 11; ++i) {
        printf("%f", pow(i, 2));
    }

    return 0;
}
```
