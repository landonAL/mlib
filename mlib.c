#ifndef MLIB_C_
#define MLIB_C_

#include <assert.h>
#include <stdbool.h>

#define PI  3.1415926535897932
#define E   2.7182818284590452
#define PHI 1.6180339887498948

double toRadiand(double deg);
float toRadianf(float deg);
double toDegreed(double deg);
float toDegreef(float deg);
int floor(double a);
int ceil(double a);
int round(double a);
double absd(double a);
float absf(float a);
double sqrt(double a);
int addi(int a, int b);
double addd(double a, double b);
float addf(float a, float b);
int subi(int a, int b);
double subd(double a, double b);
float subf(float a, float b);
int multi(int a, int b);
double multd(double a, double b);
float multf(float a, float b);
int divi(int a, int b);
double divd(double a, double b);
float divf(float a, float b);
int mod(int a, int b);
int fdiv(double a, double b);
double powd(double base, int pow);
float powf(float base, int pow);
bool isPrime(int a);

#endif // MLIB_C_

#ifdef MLIB_IMPLEMENTATION

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

int
floor(double a)
{
    return (int) a;
}

int
ceil(double a)
{
    return a > (int) a ? (int) a + 1 : (int) a;
}

int
round(double a)
{
    return a + 0.5 >= (int) a + 1
        ? (int) a + 1
        : (int) a;
}

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

double
sqrt(double a)
{
    assert(a >= 0);
    if (a == 0) return 0;

    double b = a, root;

    for (int i = 1; i < 17; ++i) {
        b = root = (b + (a / b)) / 2;
    }

    return root;
}

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

int
mod(int a, int b)
{
    assert(b > 0);
    return a % b;
}

int
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

float
powf(float base, int pow)
{
    if (pow == 0) return 1;

    float product = base;

    for (int i = 1; i < pow; ++i) product *= base;

    return product;
}

bool
isPrime(int a)
{
    if (a < 2) return false;
    if (a > 2 && a % 2 == 0) return false;

    for (int i = 2; i < a / 2; ++i) {
        if (a % i == 0) return false;
    }

    return true;
}

#endif // MLIB_IMPLEMENTATION
