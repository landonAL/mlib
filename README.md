# MLib

> [!WARNING]
> THIS LIBRARY IS A WORK IN PROGRESS! ANYTHING IS SUBJECT TO CHANGE WITHOUT NOTICE! USE THIS LIBRARY AT YOUR OWN RISK!

A mathematical library that does not have any dependencies and aims to provide advanced mathematical functionality.

The library is not concerned with displaying the image. It only fills up the memory with pixels. It's up to you what to do with those pixels.

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
