#include <assert.h>

/* VARIABLES */

float PI = 3.14159265358979324;
float E  = 2.7182818284590452;

/* ANGLES */

double
toRadiand(double deg)
{
    return deg * (PI / 180);
}

float
toRadianf(float deg)
{
    return deg * (PI / 180);
}

double
toDegreed(double rad)
{
    return rad * (180 / PI);
}

float
toDegreef(float rad)
{
    return rad * (180 / PI);
}

/* ROUNDING */

int
floor(double a)
{
    return (int) a;
}

int
ceil(double a)
{
    return (int) a + 1;
}

int
round(double a)
{
    return a + 0.5 > (int) a + 1
        ? (int) a + 1
        : (int) a;
}

/* ABSOLUTE */

double
absd(double a)
{
    return a < 0 ? a * -1 : a;
}

float
absf(float a)
{
    return a < 0 ? a * -1 : a;
}

/* SQUARE ROOT */

double
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

/* ADDITION */

int
addi(int a, int b)
{
    return a + b;
}

double
addd(double a, double b)
{
    return a + b;
}

float
addf(float a, float b)
{
    return a + b;
}

/* SUBTRACTION */

int
subi(int a, int b)
{
    return a - b;
}

double
subd(double a, double b)
{
    return a - b;
}

float
subf(float a, float b)
{
    return a - b;
}

/* MULTIPLICATION */

int
multi(int a, int b)
{
    return a * b;
}

double
multd(double a, double b)
{
    return a * b;
}

float
multf(float a, float b)
{
    return a * b;
}

/* DIVISION */

int
divi(int a, int b)
{
    assert(b > 0);
    return a / b;
}

double
divd(double a, double b)
{
    assert(b > 0);
    return a / b;
}

float
divf(float a, float b)
{
    assert(b > 0);
    return a / b;
}

/* MODULO */

int
mod(int a, int b)
{
    assert(b > 0);
    return a % b;
}

/* FLOOR DIVISION */

double
fdivd(double a, double b)
{
    assert(b > 0);
    return floord(a / b);
}

float
fdivf(float a, float b)
{
    assert(b > 0);
    return floorf(a / b);
}

/* EXPONENTS */

double
powd(double base, int pow)
{
    double product = base;

    for (int i = 1; i < pow; ++i) product *= base;

    return product;
}

float
powf(float base, int pow)
{
    float product = base;

    for (int i = 1; i < pow; ++i) product *= base;

    return product;
}
