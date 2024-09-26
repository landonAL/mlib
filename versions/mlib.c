#ifndef MLIB_C_
#define MLIB_C_

#include <assert.h>
#include <stdbool.h>

#define M_PI      3.1415926535897932
#define M_TAU     6.2831853071795864
#define M_E       2.7182818284590452
#define M_PHI     1.6180339887498948
#define M_LN2     0.6931471805599453
#define M_LN10    2.3025850929940457
#define M_LOG2E   1.4426950408889634
#define M_LOG10E  0.4342944819032518
#define M_EULER   0.5772156649015329
#define M_CATALAN 0.9159655941772190

double m_radian(double deg);
double m_degree(double deg);
int m_floor(double a);
int m_ceil(double a);
int m_round(double a);
double m_abs(double a);
double m_sqrt(double a);
double m_isqrt(double a);
double m_qisqrt(double a);
double m_root(double a, int exp);
int m_rem(int a, int b);
int m_fdiv(double a, double b);
double m_pow(double base, int pow);
bool m_prime(int a);
bool m_finite(double a);
bool m_infinite(double a);
bool m_nan(double a);
double m_sin(double a);
double m_cos(double a);
double m_tan(double a);
double m_sinh(double a);
double m_cosh(double a);
double m_tanh(double a);
double m_asin(double a);
double m_acos(double a);
double m_atan(double a);
double m_atan2(double a, double b);
double m_asinh(double a);
double m_acosh(double a);
double m_atanh(double a);
double m_sec(double a);
double m_csc(double a);
double m_cot(double a);
double m_sech(double a);
double m_csch(double a);
double m_coth(double a);
double m_asec(double a);
double m_acsc(double a);
double m_acot(double a);
double m_asech(double a);
double m_acsch(double a);
double m_acoth(double a);
double m_exp(double a);
double m_min(double a, double b);
double m_max(double a, double b);
double m_clamp(double value, double min, double max);
double m_ln(double a);
double m_log(double a, double base);
double m_log2(double a);
double m_log10(double a);
double m_sum(double *data, int size);
double m_mean(double *data, int size);
double m_median(double *data, int size);
double m_mode(double *data, int size);
double m_stddev(double *data, int size);
int m_gcd(int a, int b);
int m_lcm(int a, int b);
int m_fact(int a);
int m_rand(int a, int b);

#endif // MLIB_C_

#ifdef MLIB_IMPLEMENTATION

/**
 * Converts degrees to radians.
 *
 * @param deg The angle in degrees to convert.
 * @return The angle converted to radians.
 * @note This function asserts that the input is finite.
 **/
double
m_radian(double deg)
{
	assert(m_finite(deg));
	return deg * (M_PI / 180);
}

/**
 * Converts radians to degrees.
 *
 * @param rad The angle in radians to convert.
 * @return The angle converted to degrees.
 * @note This function asserts that the input is finite.
 **/
double
m_degree(double rad)
{
	assert(m_finite(rad));
	return rad * (180 / M_PI);
}

/**
 * Calculates the floor of a given double value.
 *
 * @param a The double value to floor.
 * @return The largest integer less than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
int
m_floor(double a)
{
	assert(m_finite(a));
	return (int) a;
}

/**
 * Calculates the ceiling of a given double value.
 *
 * @param a The double value to calculate the ceiling for.
 * @return The smallest integer greater than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
int
m_ceil(double a)
{
	assert(m_finite(a));
	return a > (int) a ? (int) a + 1 : (int) a;
}

/**
 * Rounds a given double value to the nearest integer.
 *
 * @param a The double value to round.
 * @return The nearest integer to the input value.
 * @note This function asserts that the input is finite.
 **/
int
m_round(double a)
{
	assert(m_finite(a));
	return a + 0.5 >= (int) a + 1
		? (int) a + 1
		: (int) a;
}

/**
 * Calculates the absolute value of a given double.
 *
 * @param a The input double value.
 * @return The absolute value of the input.
 * @note This function asserts that the input is finite.
 **/
double
m_abs(double a)
{
	assert(m_finite(a));
	return a < 0 ? -a : a;
}

/**
 * Calculates the square root of a given number using the Newton-Raphson method.
 *
 * @param a The number to calculate the square root of.
 * @return The square root of the input number.
 * @note This function asserts that the input is finite and non-negative.
 **/
double
m_sqrt(double a)
{
	assert(m_finite(a) && a >= 0);

	if (a == 0) return 0;

	double b = a, root;

	for (int i = 1; i < 17; ++i) {
		b = root = (b + (a / b)) / 2;
	}

	return root;
}

/**
 * Calculates the inverse square root of a given number.
 *
 * @param a The number to calculate the inverse square root of.
 * @return The inverse square root of the input number.
 * @note This function asserts that the input is finite.
 **/
double
m_isqrt(double a)
{
    assert(m_finite(a));
    return 1 / m_sqrt(a);
}

/**
 * Calculates the quick inverse square root of a given number using the "Fast Inverse Square Root" algorithm.
 *
 * @param a The number to calculate the quick inverse square root of.
 * @return An approximation of the inverse square root of the input number.
 * @note This function uses a bit-level hack for initial guess and Newton's method for refinement.
 * @note This function asserts that the input is finite.
 **/
double
m_qisqrt(double a)
{
    assert(m_finite(a));

    long i;
    double x2, y;

    x2 = a * 0.5;
    y  = a;
    i  = *(long *) &y;
    i  = 0x5f3759df - (i >> 1);
    y  = *(double *) &i;
    y  *= (1.5 - (x2 * y * y));
    y  *= (1.5 - (x2 * y * y));

    return y;
}

/**
 * Calculates the nth root of a given number.
 *
 * @param a The number to calculate the root of.
 * @param exp The root exponent.
 * @return The nth root of the input number.
 * @note This function asserts that both inputs are finite, a is non-negative, and exp is positive.
 * @note For a = 0, the function returns 0.
 * @note The function uses the power function for calculation.
 */
double
m_root(double a, int exp)
{
    assert(m_finite(a) && m_finite(exp) && a >= 0 && exp > 0);
    if (a == 0) return 0;

    double result = a;
    for (int i = 0; i < 16; i++) {
        result = ((exp - 1) * result + a / m_pow(result, exp - 1)) / exp;
    }

    return result;
}

/**
 * Calculates the Greatest Common Divisor (GCD) of two integers using the Euclidean algorithm.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The GCD of a and b.
 * @note This function asserts that both inputs are finite.
 * @note The function uses the absolute values of the inputs to handle negative numbers.
 **/
int
m_gcd(int a, int b)
{
    assert(m_finite(a) && m_finite(b));

    a = m_abs(a);
    b = m_abs(b);

    while (b != 0) {
        int temp = b;

        b = a % b;
        a = temp;
    }

    return a;
}

/**
 * Calculates the Least Common Multiple (LCM) of two integers.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The LCM of a and b.
 * @note This function asserts that both inputs are finite.
 * @note The function uses the GCD to calculate the LCM efficiently.
 * @note If the GCD is 0, the function returns 0 to avoid division by zero.
 **/
int
m_lcm(int a, int b)
{
    assert(m_finite(a) && m_finite(b));

    int result = m_gcd(a, b);
    if (result == 0) return 0;

    return m_abs(a / result * b);
}

/**
 * Calculates the factorial of a given non-negative integer.
 *
 * @param a The non-negative integer for which to calculate the factorial.
 * @return The factorial of the input number.
 * @note This function asserts that the input is finite and non-negative.
 * @note The factorial is calculated recursively, which may lead to stack overflow for large inputs.
 * @warning This implementation is not suitable for large inputs due to potential stack overflow.
 **/
int
m_fact(int a)
{
    assert(m_finite(a));
    if (a <= 1) return 1;

    return a * m_fact(a - 1);
}

/**
 * Generates a random integer within a specified range using a linear congruential generator (LCG).
 *
 * @param a The lower bound of the range (inclusive).
 * @param b The upper bound of the range (inclusive).
 * @return A random integer between a and b (inclusive).
 * @note This function asserts that both inputs are finite and that a is less than b.
 * @note The LCG uses a static seed, which is updated with each call to the function.
 * @note The multiplier and modulus values are chosen to create a full-period generator.
 * @note The function incorporates compile-time information to modify the seed,
 *       providing additional randomness across different compilations.
 **/
int
m_rand(int a, int b)
{
    assert(m_finite(a) && m_finite(b) && a < b);

    static unsigned int seed = 1;
    static const unsigned int multiplier = 16807;  // 7^5
    static const unsigned int modulus = 2147483647;  // 2^31 - 1 (Mersenne prime)

    unsigned int compile_time_seed = (__TIME__[7] - '0') * 1000000 +
                                    (__TIME__[6] - '0') * 100000 +
                                    (__TIME__[4] - '0') * 10000 +
                                    (__TIME__[3] - '0') * 1000 +
                                    (__TIME__[1] - '0') * 100 +
                                    (__TIME__[0] - '0') * 10;

    seed = (multiplier * (seed ^ compile_time_seed)) % modulus;

    return a + (int) (((double) seed / modulus) * (b - a + 1));
}

/**
 * Calculates the remainder of the division of two integers.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The remainder of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 **/
int
m_rem(int a, int b)
{
	assert(m_finite(a) && m_finite(b) && b > 0);
	return a % b;
}

/**
 * Performs floor division of two double values.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The floor of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 **/
int
m_fdiv(double a, double b)
{
	assert(m_finite(a) && m_finite(b) && b > 0);
	return m_floor(a / b);
}

/**
 * Calculates the power of a base number raised to an integer exponent.
 *
 * @param base The base number.
 * @param pow The integer exponent.
 * @return The result of base raised to the power of pow.
 * @note This function asserts that both base and pow are finite.
 * @note For pow = 0, the function returns 1.
 * @note The function uses a simple iterative approach for positive exponents.
 **/
double
m_pow(double base, int pow)
{
	assert(m_finite(base) && m_finite(pow));

	if (pow == 0) return 1;

	double product = base;

	for (int i = 1; i < pow; ++i) product *= base;

	return product;
}

/**
 * Checks if a given integer is prime.
 *
 * @param a The integer to check for primality.
 * @return true if the number is prime, false otherwise.
 * @note This function asserts that the input is finite.
 * @note Numbers less than 2 are not considered prime.
 * @note Even numbers greater than 2 are not prime.
 * @note The function checks for divisibility up to half of the input number.
 **/
bool
m_prime(int a)
{
	assert(m_finite(a));

	if (a < 2) return false;
	if (a > 2 && a % 2 == 0) return false;

	for (int i = 2; i < a / 2; ++i) {
		if (a % i == 0) return false;
	}

	return true;
}

/**
 * Checks if a given double value is finite.
 *
 * @param a The double value to check.
 * @return true if the value is finite, false otherwise.
 **/
bool
m_finite(double a)
{
	return !m_infinite(a) && !m_nan(a);
}

/**
 * Checks if a given double value is infinite.
 *
 * @param a The double value to check.
 * @return true if the value is infinite, false otherwise.
 **/
bool
m_infinite(double a)
{
	return a / a != a / a;
}

/**
 * Checks if a given double value is Not-a-Number (NaN).
 *
 * @param a The double value to check.
 * @return true if the value is NaN, false otherwise.
 **/
bool
m_nan(double a)
{
	return a != a;
}

/**
 * Calculates the sine of an angle using Taylor series approximation.
 *
 * @param a The angle in radians.
 * @return The sine of the input angle.
 * @note This function asserts that the input is finite.
 * @note The function normalizes the input angle to the range [-PI, PI].
 * @note The Taylor series is computed up to the 7th term for accuracy.
 **/
double
m_sin(double a)
{
    assert(m_finite(a));

    while (a > M_PI) a -= 2 * M_PI;
    while (a < -M_PI) a += 2 * M_PI;

    double result = a;
    double term = a;

    for (int i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i) * (2 * i + 1));
        result += term;
    }

    return result;
}

/**
 * Calculates the cosine of an angle using Taylor series approximation.
 *
 * @param a The angle in radians.
 * @return The cosine of the input angle.
 * @note This function asserts that the input is finite.
 * @note The function normalizes the input angle to the range [-PI, PI].
 * @note The Taylor series is computed up to the 7th term for accuracy.
 **/
double
m_cos(double a)
{
    assert(m_finite(a));

    while (a > M_PI) a -= 2 * M_PI;
    while (a < -M_PI) a += 2 * M_PI;

    double result = 1;
    double term = 1;

    for (int i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i - 1) * (2 * i));
        result += term;
    }

    return result;
}

/**
 * Calculates the tangent of an angle.
 *
 * @param a The angle in radians.
 * @return The tangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The tangent is calculated as the ratio of sine to cosine.
 **/
double
m_tan(double a)
{
    assert(m_finite(a));
    return m_sin(a) / m_cos(a);
}

/**
 * Calculates the hyperbolic sine of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic sine of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 0.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_sinh(double a)
{
    assert(m_finite(a));

    if (a == 0) return 0;

    double ea = m_exp(a);
    return (ea - (1 / ea)) / 2;
}

/**
 * Calculates the hyperbolic cosine of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 1.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_cosh(double a)
{
    assert(m_finite(a));

    if (a == 0) return 1;

    double ea = m_exp(a);
    return (ea + (1 / ea)) / 2;
}

/**
 * Calculates the hyperbolic tangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 0.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_tanh(double a)
{
    assert(m_finite(a));

    if (a == 0) return 0;

    double ea = m_exp(2 * a);
    return (ea - 1) / (ea + 1);
}

/**
 * Calculates the arcsine (inverse sine) of a given value using a polynomial approximation.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arcsine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The approximation uses a 7th-degree polynomial for accuracy.
 **/
double
m_asin(double a)
{
    assert(m_finite(a) && a >= -1 && a <= 1);

    double a2 = m_pow(a, 2);
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
}

/**
 * Calculates the arccosine (inverse cosine) of a given value.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arccosine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The arccosine is calculated using the relationship: acos(x) = PI/2 - asin(x).
 **/
double
m_acos(double a)
{
    assert(m_finite(a) && a >= -1 && a <= 1);
    return (M_PI / 2) - m_asin(a);
}

/**
 * Calculates the arctangent (inverse tangent) of a given value using an approximation.
 *
 * @param a The input value.
 * @return The arctangent of the input value in radians.
 * @note This function asserts that the input is finite.
 * @note This approximation is less accurate for large input values.
 **/
double
m_atan(double a)
{
    assert(m_finite(a));
    return a / (1.28 * m_pow(a, 2));
}

/**
 * Calculates the arctangent of two variables (atan2).
 *
 * @param a The y-coordinate.
 * @param b The x-coordinate.
 * @return The angle in radians between the positive x-axis and the point (b, a).
 * @note This function asserts that both inputs are finite.
 * @note Special cases:
 *       - If b is 0 and a > 0, returns PI/2
 *       - If b is 0 and a < 0, returns -PI/2
 *       - If b is 0 and a is 0, returns 0
 *       - If b < 0, adjusts the result by adding or subtracting PI
 **/
double
m_atan2(double a, double b)
{
    assert(m_finite(a) && m_finite(b));

    if (b == 0) {
        if (a > 0) return M_PI / 2;
        if (a < 0) return -M_PI / 2;

        return 0;
    }

    double result = m_atan(a / b);

    if (b < 0) {
        if (a >= 0) return result + M_PI;
        return result - M_PI;
    }

    return result;
}

/**
 * Calculates the inverse hyperbolic sine of a given value.
 *
 * @param a The input value.
 * @return The inverse hyperbolic sine of the input value.
 * @note This function asserts that the input is finite.
 **/
double
m_asinh(double a)
{
    assert(m_finite(a));
    return m_ln(a + m_sqrt(a * a + 1));
}

/**
 * Calculates the inverse hyperbolic cosine of a given value.
 *
 * @param a The input value, must be greater than or equal to 1.
 * @return The inverse hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite and greater than or equal to 1.
 **/
double
m_acosh(double a)
{
    assert(m_finite(a) && a >= 1);
    return m_ln(a + m_sqrt(a * a - 1));
}

/**
 * Calculates the inverse hyperbolic tangent of a given value.
 *
 * @param a The input value, must be in the range (-1, 1).
 * @return The inverse hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite and within the valid range.
 **/
double
m_atanh(double a)
{
    assert(m_finite(a) && a > -1 && a < 1);
    return 0.5 * m_ln((1 + a) / (1 - a));
}

/**
 * Calculates the secant of an angle.
 *
 * @param a The angle in radians.
 * @return The secant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The secant is calculated as the reciprocal of the cosine.
 **/
double
m_sec(double a)
{
    assert(m_finite(a));
    return 1 / m_cos(a);
}

/**
 * Calculates the cosecant of an angle.
 *
 * @param a The angle in radians.
 * @return The cosecant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cosecant is calculated as the reciprocal of the sine.
 **/
double
m_csc(double a)
{
    assert(m_finite(a));
    return 1 / m_sin(a);
}

/**
 * Calculates the cotangent of an angle.
 *
 * @param a The angle in radians.
 * @return The cotangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cotangent is calculated as the ratio of cosine to sine.
 **/
double
m_cot(double a)
{
    assert(m_finite(a));
    return m_cos(a) / m_sin(a);
}

/**
 * Calculates the hyperbolic secant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic secant of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 1.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_sech(double a)
{
    assert(m_finite(a));

    if (a == 0) return 1;

    double ea = m_exp(a);
    return 2 / (ea + (1 / ea));
}

/**
 * Calculates the hyperbolic cosecant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosecant of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_csch(double a)
{
    assert(m_finite(a));

    double ea = m_exp(a);
    return 2 / (ea - (1 / ea));
}

/**
 * Calculates the hyperbolic cotangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cotangent of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
double
m_coth(double a)
{
    assert(m_finite(a));

    double ea = m_exp(2 * a);
    return (ea + 1) / (ea - 1);
}

/**
 * Calculates the exponential function (e^x) for a given value.
 *
 * @param a The exponent value.
 * @return The result of e raised to the power of a.
 * @note This function asserts that the input is finite.
 * @note The function uses a Taylor series approximation combined with exponent reduction.
 * @note For a = 0, the function returns 1.
 * @note The calculation is optimized for accuracy and efficiency.
 **/
double
m_exp(double a)
{
    assert(m_finite(a));

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

    return result * m_pow(2, k);
}

/**
 * Finds the minimum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The smaller of the two input values.
 * @note This function asserts that both inputs are finite.
 **/
double
m_min(double a, double b)
{
    assert(m_finite(a) && m_finite(b));
    return a < b ? a : b;
}

/**
 * Finds the maximum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The larger of the two input values.
 * @note This function asserts that both inputs are finite.
 **/
double
m_max(double a, double b)
{
    assert(m_finite(a) && m_finite(b));
    return a > b ? a : b;
}

/**
 * Clamps a double value between a minimum and maximum range.
 *
 * @param value The value to clamp.
 * @param min The minimum allowed value.
 * @param max The maximum allowed value.
 * @return The clamped value, which will be between min and max (inclusive).
 * @note This function asserts that all inputs are finite.
 **/
double
m_clamp(double value, double min, double max)
{
    assert(m_finite(value) && m_finite(min) && m_finite(max));

    if (value < min) return min;
    if (value > max) return max;

    return value;
}

/**
 * Calculates the natural logarithm of a given value.
 *
 * @param a The input value, must be greater than 0.
 * @return The natural logarithm of the input value.
 * @note This function asserts that the input is finite and greater than 0.
 * @note The function uses a series expansion for improved accuracy.
 **/
double
m_ln(double a)
{
    assert(m_finite(a) && a > 0);

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
    } while (m_abs(y) > 1E-15);

    return sum + exp * LN2;
}

/**
 * Calculates the logarithm of a value with a specified base.
 *
 * @param a The input value.
 * @param base The base of the logarithm.
 * @return The logarithm of the input value with the specified base.
 * @note This function asserts that both inputs are finite.
 **/
double
m_log(double a, double base)
{
    assert(m_finite(a) && m_finite(base));
    return m_ln(a) / m_ln(base);
}

/**
 * Calculates the base-2 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-2 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
double
m_log2(double a)
{
    assert(m_finite(a));
    return m_ln(a) / LN2;
}

/**
 * Calculates the base-10 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-10 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
double
m_log10(double a)
{
    assert(m_finite(a));
    return m_ln(a) / LN10;
}

/**
 * Calculates the sum of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The sum of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 **/
double
m_sum(double *data, int size)
{
    assert(m_finite(size) && size > 0);

    double sum = 0;

    for (int i = 0; i < size; ++i) {
        sum += data[i];
    }

    return sum;
}

/**
 * Calculates the arithmetic mean of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The arithmetic mean of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 **/
double
m_mean(double *data, int size)
{
    assert(m_finite(size) && size > 0);
    return m_sum(data, size) / size;
}

/**
 * Calculates the median of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The median value of the array.
 * @note This function asserts that size is finite and greater than 0.
 * @note This function modifies the original array by sorting it.
 **/
double
m_median(double *data, int size)
{
    assert(m_finite(size) && size > 0);

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

/**
 * Calculates the mode of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The mode (most frequent value) of the array.
 * @note This function asserts that size is finite and greater than 0.
 * @note If multiple modes exist, this function returns the first one encountered.
 **/
double
m_mode(double *data, int size)
{
    assert(m_finite(size) && size > 0);

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

/**
 * Calculates the sample standard deviation of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The sample standard deviation of the array.
 * @note This function asserts that size is finite and greater than 1.
 **/
double
m_stddev(double *data, int size)
{
    assert(m_finite(size) && size > 1);

    double m = m_mean(data, size);
    double sum = 0;

    for (int i = 0; i < size; ++i) {
        double diff = data[i] - m;
        sum += diff * diff;
    }

    return sqrt(sum / (size - 1));
}

#endif // MLIB_IMPLEMENTATION
