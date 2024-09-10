#ifndef MLIB_C_
#define MLIB_C_

#include <assert.h>
#include <stdbool.h>

#define PI  3.1415926535897932
#define E   2.7182818284590452
#define PHI 1.6180339887498948

typedef struct
{
    double x;
    double y;
} Vector2;

typedef struct
{
    double x;
    double y;
    double z;
} Vector3;

double toRadiand(double deg);
double toDegreed(double deg);
int floor(double a);
int ceil(double a);
int round(double a);
double abs(double a);
double sqrt(double a);
int mod(int a, int b);
int fdiv(double a, double b);
double pow(double base, int pow);
bool isPrime(int a);
void v2add(Vector2 a, Vector2 b);
void v2sub(Vector2 a, Vector2 b);
void v2mult(Vector2 a, Vector2 b);
void v2div(Vector2 a, Vector2 b);
void v2perp(Vector2 a);
double v2cross(Vector2 a, Vector2 b);
double v2project(Vector2 a);
void v2reverse(Vector2 a);
void v2set(Vector2 a, Vector2 b);
void v2zero(Vector2 a);
double v2length(Vector2 a);
double v2lengthSq(Vector2 a);
double v2distance(Vector2 a, Vector2 b);
double v2distanceSq(Vector2 a, Vector2 b);
void v2norm(Vector2 a);
void v3add(Vector3 a, Vector3 b);
void v3sub(Vector3 a, Vector3 b);
void v3mult(Vector3 a, Vector3 b);
void v3div(Vector3 a, Vector3 b);
double v3cross(Vector3 a, Vector3 b);
double v3project(Vector3 a);
void v3reverse(Vector3 a);
void v3set(Vector3 a, Vector3 b);
void v3zero(Vector3 a);
double v3length(Vector3 a);
double v3lengthSq(Vector3 a);
double v3distance(Vector3 a, Vector3 b);
double v3distanceSq(Vector3 a, Vector3 b);
void v3norm(Vector3 a);

#endif // MLIB_C_

#ifdef MLIB_IMPLEMENTATION

double
toRadiand(double deg)
{
    return deg * (PI / 180);
}

double
toDegreed(double rad)
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
abs(double a)
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
pow(double base, int pow)
{
    if (pow == 0) return 1;

    double product = base;

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

void
v2add(Vector2 a, Vector2 b)
{
    a.x += b.x;
    a.y += b.y;
}

void
v2sub(Vector2 a, Vector2 b)
{
    a.x -= b.x;
    a.y -= b.y;
}

void
v2mult(Vector2 a, Vector2 b)
{
    a.x *= b.x;
    a.y *= b.y;
}

void
v2div(Vector2 a, Vector2 b)
{
    a.x /= b.x;
    a.y /= b.y;
}

void
v2perp(Vector2 a)
{
    a.x = -a.y;
    a.y = a.x;
}

double
v2cross(Vector2 a, Vector2 b)
{
    return (a.x * b.y) - (a.y * b.x);
}

double
v2project(Vector2 a)
{
    return v2lengthSq(a) / v2length(a);
}

void
v2reverse(Vector2 a)
{
    a.x = -a.x;
    a.y = -a.y;
}

void
v2set(Vector2 a, Vector2 b)
{
    a.x = b.x;
    a.y = b.y;
}

void
v2zero(Vector2 a)
{
    a.x = a.y = 0;
}

double
v2length(Vector2 a)
{
    return sqrt((a.x * a.x) + (a.y * a.y));
}

double
v2lengthSq(Vector2 a)
{
    return (a.x * a.x) + (a.y * a.y);
}

double
v2distance(Vector2 a, Vector2 b)
{
    a.x -= b.x;
    a.y -= b.y;

    return sqrt((a.x * a.x) + (a.y * a.y));
}

double
v2distanceSq(Vector2 a, Vector2 b)
{
    a.x -= b.x;
    a.y -= b.y;

    return (a.x * a.x) + (a.y * a.y);
}

void
v2norm(Vector2 a)
{
    double magnitude = v2length(a);

    a.x /= magnitude;
    a.y /= magnitude;
}

void
v3add(Vector3 a, Vector3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

void
v3sub(Vector3 a, Vector3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
}

void
v3mult(Vector3 a, Vector3 b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
}

void
v3div(Vector3 a, Vector3 b)
{
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
}

double
v3cross(Vector3 a, Vector3 b)
{
    return (a.x * b.y) - (a.y * b.x);
}

double
v3project(Vector3 a)
{
    return v3lengthSq(a) / v3length(a);
}

void
v3reverse(Vector3 a)
{
    a.x = -a.x;
    a.y = -a.y;
    a.z = -a.z;
}

void
v3set(Vector3 a, Vector3 b)
{
    a.x = b.x;
    a.y = b.y;
    a.z = b.z;
}

void
v3zero(Vector3 a)
{
    a.x = a.y = a.z = 0;
}

double
v3length(Vector3 a)
{
    return sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}

double
v3lengthSq(Vector3 a)
{
    return (a.x * a.x) + (a.y * a.y) + (a.z * a.z);
}

double
v3distance(Vector3 a, Vector3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;

    return sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}

double
v3distanceSq(Vector3 a, Vector3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;

    return (a.x * a.x) + (a.y * a.y) + (a.z * a.z);
}

void
v3norm(Vector3 a)
{
    double magnitude = v3length(a);

    a.x /= magnitude;
    a.y /= magnitude;
    a.z /= magnitude;
}

#endif // MLIB_IMPLEMENTATION
