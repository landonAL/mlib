#ifndef MLIB_C_
#define MLIB_C_

#include <assert.h>

#ifndef MLIBDEF
#define MLIBDEF static inline
#endif

#define PI 3.14159265358979324
#define E  2.7182818284590452

MLIBDEF double toRadiand(double deg);
MLIBDEF float toRadianf(float deg);
MLIBDEF double toDegreed(double deg);
MLIBDEF float toDegreef(float deg);
MLIBDEF int floor(double a);
MLIBDEF int ceil(double a);
MLIBDEF int round(double a);
MLIBDEF double absd(double a);
MLIBDEF float absf(float a);
MLIBDEF double sqrt(double a);
MLIBDEF int addi(int a, int b);
MLIBDEF double addd(double a, double b);
MLIBDEF float addf(float a, float b);
MLIBDEF int subi(int a, int b);
MLIBDEF double subd(double a, double b);
MLIBDEF float subf(float a, float b);
MLIBDEF int multi(int a, int b);
MLIBDEF double multd(double a, double b);
MLIBDEF float multf(float a, float b);
MLIBDEF int divi(int a, int b);
MLIBDEF double divd(double a, double b);
MLIBDEF float divf(float a, float b);
MLIBDEF int mod(int a, int b);
MLIBDEF int fdiv(double a, double b);
MLIBDEF double powd(double base, int pow);
MLIBDEF float powf(float base, int pow);

#endif // MLIB_C_

#ifdef MLIB_IMPLEMENTATION

MLIBDEF double
toRadiand(double deg)
{
    return deg * (PI / 180);
}

MLIBDEF float
toRadianf(float deg)
{
    return deg * (PI / 180);
}

MLIBDEF double
toDegreed(double rad)
{
    return rad * (180 / PI);
}

MLIBDEF float
toDegreef(float rad)
{
    return rad * (180 / PI);
}

MLIBDEF int
floor(double a)
{
    return (int) a;
}

MLIBDEF int
ceil(double a)
{
    return a > (int) a ? (int) a + 1 : (int) a;
}

MLIBDEF int
round(double a)
{
    return a + 0.5 >= (int) a + 1
        ? (int) a + 1
        : (int) a;
}

MLIBDEF double
absd(double a)
{
    return a < 0 ? a * -1 : a;
}

MLIBDEF float
absf(float a)
{
    return a < 0 ? a * -1 : a;
}

MLIBDEF double
sqrt(double a)
{
    double b = a, root;

    for (int i = 0; i < 5; ++i) {
        root = 0.5 * (b + (a / b));

	if ((root - b) < 1)
	    break;

	b = root;
    }

    return root;
}

MLIBDEF int
addi(int a, int b)
{
    return a + b;
}

MLIBDEF double
addd(double a, double b)
{
    return a + b;
}

MLIBDEF float
addf(float a, float b)
{
    return a + b;
}

MLIBDEF int
subi(int a, int b)
{
    return a - b;
}

MLIBDEF double
subd(double a, double b)
{
    return a - b;
}

MLIBDEF float
subf(float a, float b)
{
    return a - b;
}

MLIBDEF int
multi(int a, int b)
{
    return a * b;
}

MLIBDEF double
multd(double a, double b)
{
    return a * b;
}

MLIBDEF float
multf(float a, float b)
{
    return a * b;
}

MLIBDEF int
divi(int a, int b)
{
    assert(b > 0);
    return a / b;
}

MLIBDEF double
divd(double a, double b)
{
    assert(b > 0);
    return a / b;
}

MLIBDEF float
divf(float a, float b)
{
    assert(b > 0);
    return a / b;
}

MLIBDEF int
mod(int a, int b)
{
    assert(b > 0);
    return a % b;
}

MLIBDEF int
fdiv(double a, double b)
{
    assert(b > 0);
    return floor(a / b);
}

double
powd(double base, int pow)
{
    if (pow == 0) return 1;

    double product = base;

    for (int i = 1; i < pow; ++i) product *= base;

    return product;
}

MLIBDEF float
powf(float base, int pow)
{
    if (pow == 0) return 1;

    float product = base;

    for (int i = 1; i < pow; ++i) product *= base;

    return product;
}

#endif // MLIB_IMPLEMENTATION
