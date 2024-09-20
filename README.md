# MLib

> [!WARNING]
> This library is a work in progress, and anything is subject to change without notice. Use this library at your own risk!

A mathematical library that does not have any dependencies and aims to provide advanced mathematical functionality.

The library itself does not need to be built. You can simply copy-paste mlib into your project directory and include it.

In its primary language, MLib is a classic [stb-style](https://github.com/nothings/stb) single header library, meaning that by default it acts like a header, but if you add `#define MLIB_IMPLEMENTATION` prior to including the library it acts like a source file.

## Example
### C
```c
#include <stdio.h>

#define MLIB_IMPLEMENTATION
#include "mlib.c"

int
main()
{
    for (int i = 1; i < 11; ++i) {
        printf("%f\n", pow(i, 2));
    }

    return 0;
}
```

### Python
```python
import mlib

def main():
    for i in range(11):
        print(pow(i, 2))

if __name__ == "__main__":
    main()
```

## Known issues
None currently found
    
