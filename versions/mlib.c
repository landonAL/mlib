#ifndef MLIB_C_
#define MLIB_C_

#include <assert.h>
#include <stdbool.h>

#define PI     3.1415926535897932
#define TAU    6.2831853071795864
#define E      2.7182818284590452
#define PHI    1.6180339887498948
#define LN2    0.6931471805599453
#define LN10   2.3025850929940457
#define LOG2E  1.4426950408889634
#define LOG10E 0.4342944819032518
#define EULER  0.5772156649015329

double toRadian(double deg);
double toDegree(double deg);
int floor(double a);
int ceil(double a);
int round(double a);
double abs(double a);
double sqrt(double a);
double isqrt(double a);
double qisqrt(double a);
int mod(int a, int b);
int fdiv(double a, double b);
double pow(double base, int pow);
bool isPrime(int a);
bool isFinite(double a);
bool isInfinite(double a);
bool isNaN(double a);
double sin(double a);
double cos(double a);
double tan(double a);
double sinh(double a);
double cosh(double a);
double tanh(double a);
double asin(double a);
double acos(double a);
double atan(double a);
double atan2(double a, double b);
double asinh(double a);
double acosh(double a);
double atanh(double a);
double sec(double a);
double csc(double a);
double cot(double a);
double exp(double a);
double min(double a, double b);
double max(double a, double b);
double clamp(double value, double min, double max);
double ln(double a);
double log2(double a);
double log10(double a);
double sum(double *data, int size);
double mean(double *data, int size);
double median(double *data, int size);
double mode(double *data, int size);
double stddev(double *data, int size);
int gcd(int a, int b);
int lcm(int a, int b);
int fact(int a);

#endif // MLIB_C_

#ifdef MLIB_IMPLEMENTATION

double
toRadian(double deg)
{
	return deg * (PI / 180);
}

double
toDegree(double rad)
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
	return a < 0 ? -a : a;
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

double
isqrt(double a)
{
    return 1 / sqrt(a);
}

double
qisqrt(double a)
{
    long i;
    float x2, y;

    x2 = a * 0.5;
    y  = a;
    i  = *(long *) &y;
    i  = 0x5f3759df - (i >> 1);
    y  = *(float *) &i;
    y  *= (1.5 - (x2 * y * y));
    y  *= (1.5 - (x2 * y * y));

    return y;
}

int
gcd(int a, int b)
{
    a = abs(a);
    b = abs(b);

    while (b != 0) {
        int temp = b;

        b = a % b;
        a = temp;
    }

    return a;
}

int
lcm(int a, int b)
{
    int result = gcd(a, b);
    if (result == 0) return 0;

    return abs(a / result * b);
}

int
fact(int a)
{
    int result = 1;

    for (int i = 2; i <= a; ++i) {
        result *= i;
    }

    return result;
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

bool
isFinite(double a)
{
	return !isInfinite(a) && !isNaN(a);
}

bool
isInfinite(double a)
{
	return a / a != a / a;
}

bool
isNaN(double a)
{
	return a != a;
}

double
sin(double a)
{
    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;

    double result = a;
    double term = a;

    for (int i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i) * (2 * i + 1));
        result += term;
    }

    return result;
}

double
cos(double a)
{
    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;

    double result = 1;
    double term = 1;

    for (int i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i - 1) * (2 * i));
        result += term;
    }

    return result;
}

double
tan(double a)
{
    double s = sin(a);
    double c = cos(a);

    if (c == 0) return isInfinite(a);

    return s / c;
}

double
sinh(double a)
{
    if (a == 0) return 0;
    if (isInfinite(a)) return a;

    double ea = exp(a);
    return (ea - (1 / ea)) / 2;
}

double
cosh(double a)
{
    if (a == 0) return 1;
    if (isInfinite(a)) return abs(a);

    double ea = exp(a);
    return (ea + (1 / ea)) / 2;
}

double
tanh(double a)
{
    if (a == 0) return 0;
    if (isInfinite(a)) return a > 0 ? 1 : -1;

    double ea = exp(2 * a);
    return (ea - 1) / (ea + 1);
}

double
asin(double a)
{
    assert(a >= -1 && a <= 1);

    double a2 = pow(a, 2);
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
}

double
acos(double a)
{
    assert(a >= -1 && a <= 1);
    return (PI / 2) - asin(a);
}

double
atan(double a)
{
    return a / (1.28 * pow(a, 2));
}

double
atan2(double a, double b)
{
    if (b == 0) {
        if (a > 0) return PI / 2;
        if (a < 0) return -PI / 2;

        return 0;
    }

    double result = atan(a / b);

    if (b < 0) {
        if (a >= 0) return result + PI;
        return result - PI;
    }

    return result;
}

double
asinh(double a)
{
    return ln(a + sqrt(a * a + 1));
}

double
acosh(double a)
{
    assert(a >= 1);
    return ln(a + sqrt(a * a - 1));
}

double
atanh(double a)
{
    assert(a > -1 && a < 1);
    return 0.5 * ln((1 + a) / (1 - a));
}

double
sec(double a)
{
    double c = cos(a);
    if (c == 0) return isInfinite(a);

    return 1 / c;
}

double
csc(double a)
{
    double s = sin(a);
    if (s == 0) return isInfinite(a);

    return 1 / s;
}

double
cot(double a)
{
    double s = sin(a);
    double c = cos(a);

    if (s == 0) return isInfinite(a);

    return c / s;
}

double
exp(double a)
{
    assert(isFinite(a));
    if (a == 0) return 1;

    int k = (int) (a * LOG2E);
    double r = a - k * LN2;
    double result = r + 1;
    double term = r;

    for (int i = 2; i <= 12; ++i) {
        term *= r / i;
        result += term;

        if (term < 1E-15 * result) break;
    }

    return result * pow(2, k);
}

double
min(double a, double b)
{
    return a < b ? a : b;
}

double
max(double a, double b)
{
    return a > b ? a : b;
}

double
clamp(double value, double min, double max)
{
    if (value < min) return min;
    if (value > max) return max;

    return value;
}

double
ln(double a)
{
    assert(a > 0);
    if (a == 1) return 0;

    int exp = 0;

    while (a > 2) { a /= 2; exp++; }
    while (a < 1) { a *= 2; exp--; }

    a--;

    double y = a;
    double sum = y;
    int i = 1;

    do {
        i++;
        y *= -a * (i - 1) / i;
        sum += y;
    } while (abs(y) > 1E-15);

    return sum + exp * LN2;
}

double
log2(double a)
{
    return ln(a) / LN2;
}

double
log10(double a)
{
    return ln(a) / LN10;
}

double
sum(double *data, int size)
{
    assert(size > 0);

    double sum = 0;

    for (int i = 0; i < size; ++i) {
        sum += data[i];
    }

    return sum;
}

double
mean(double *data, int size)
{
    return sum(data) / size;
}

double
median(double *data, int size)
{
    assert(size > 0);

    for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size - i - 1; ++j) {
            if (data[j] > data[j + 1]) {
                double temp = data[j];

                data[j] = data[j + 1];
                data[j + 1] = temp;
            }
        }
    }

    if (size % 2 == 0) {
        return (data[(size / 2) - 1] + data[size / 2]) / 2;
    } else {
        return data[size / 2];
    }
}

double
mode(double *data, int size)
{
    assert(size > 0);

    double mode = data[0];
    int maxCount = 1, currentCount = 1;

    for (int i = 1; i < size; ++i) {
        if (data[i] == data[i - 1]) {
            ++currentCount;
        } else {
            if (currentCount > maxCount) {
                maxCount = currentCount;
                mode = data[i - 1];
            }

            currentCount = 1;
        }
    }

    if (currentCount > maxCount) {
        mode = data[size - 1];
    }

    return mode;
}

double
stddev(double *data, int size)
{
    assert(size > 1);

    double m = mean(data, size);
    double sum = 0;

    for (int i = 0; i < size; i++) {
        double diff = data[i] - m;
        sum += diff * diff;
    }

    return sqrt(sum / (size - 1));
}

#endif // MLIB_IMPLEMENTATION
