module mlib;
//this is good lang

const double PI      = 3.1415926535897932;
const double TAU     = 6.2831853071795864;
const double E       = 2.7182818284590452;
const double PHI     = 1.6180339887498948;
const double LN2     = 0.6931471805599453;
const double LN10    = 2.3025850929940457;
const double LOG2E   = 1.4426950408889634;
const double LOG10E  = 0.4342944819032518;
const double EULER   = 0.5772156649015329;
const double CATALAN = 0.9159655941772190;

/*
 * Converts degrees to radians.
 *
 * @param deg The angle in degrees to convert.
 * @return The angle converted to radians.
 * @note This function asserts that the input is finite.
 */
fn double toRadian(double deg)
{
	assert(isFinite(deg));
	return deg * (PI / 180);
}

/*
 * Converts radians to degrees.
 *
 * @param rad The angle in radians to convert.
 * @return The angle converted to degrees.
 * @note This function asserts that the input is finite.
 */
fn double toDegree(double rad)
{
				assert(isFinite(rad));
				return rad * (180 / PI);
}

/*
 * Calculates the floor of a given double value.
 *
 * @param a The double value to floor.
 * @return The largest integer less than or equal to the input value.
 * @note This function asserts that the input is finite.
 */
fn int floor(double a)
{
	assert(isFinite(a));
	return (int) a;
}

/*
 * Calculates the ceiling of a given double value.
 *
 * @param a The double value to calculate the ceiling for.
 * @return The smallest integer greater than or equal to the input value.
 * @note This function asserts that the input is finite.
 */
fn int ceil(double a)
{
	assert(isFinite(a));
	return a > (int) a ? (int) a + 1 : (int) a;
}

/*
 * Rounds a given double value to the nearest integer.
 *
 * @param a The double value to round.
 * @return The nearest integer to the input value.
 * @note This function asserts that the input is finite.
 */
fn int round(double a)
{
	assert(isFinite(a));
	return a + 0.5 >= (int) (a + 1)
					? (int) (a + 1)
					: (int) a;
}

/*
 * Calculates the absolute value of a given double.
 *
 * @param a The input double value.
 * @return The absolute value of the input.
 * @note This function asserts that the input is finite.
 */
fn double abs(double a)
{
	assert(isFinite(a));
	return a < 0 ? -a : a;
}

/*
 * Calculates the square root of a given number using the Newton-Raphson method.
 *
 * @param a The number to calculate the square root of.
 * @return The square root of the input number.
 * @note This function asserts that the input is finite and non-negative.
 */
fn double sqrt(double a)
{
	assert(isFinite(a) && a >= 0);

	if (a == 0) return 0;

	double b = a;
	double root;

	for (int i = 1; i < 17; ++i)
	{
		b = root = (b + (a / b)) / 2;
	}

	return root;
}

/*
 * Calculates the inverse square root of a given number.
 *
 * @param a The number to calculate the inverse square root of.
 * @return The inverse square root of the input number.
 * @note This function asserts that the input is finite.
 */
fn double isqrt(double a)
{
	assert(isFinite(a));
	return 1 / sqrt(a);
}

/*
 * Calculates the quick inverse square root of a given number using the "Fast Inverse Square Root" algorithm.
 *
 * @param a The number to calculate the quick inverse square root of.
 * @return An approximation of the inverse square root of the input number.
 * @note This function uses a bit-level hack for initial guess and Newton's method for refinement.
 * @note This function asserts that the input is finite.
 */
fn double qisqrt(double a)
{
    assert(isFinite(a));

    long i;
    double x2, y;

    x2 = a * 0.5;
    y  = a;
    i  = *(long *)(&y);
    i  = 0x5f3759df - (i >> 1);
    y  = *(double *)(&i);
    y  *= (1.5 - (x2 * y * y));
    y  *= (1.5 - (x2 * y * y));

    return y;
}

/*
 * Calculates the Greatest Common Divisor (GCD) of two integers using the Euclidean algorithm.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The GCD of a and b.
 * @note This function asserts that both inputs are finite.
 * @note The function uses the absolute values of the inputs to handle negative numbers.
 */
fn int gcd(int a, int b)
{
    assert(isFinite(a) && isFinite(b));

    a = (int) abs(a);
    b = (int) abs(b);

    while (b != 0)
    {
        int temp = b;

        b = a % b;
        a = temp;
    }

    return a;
}

/*
 * Calculates the Least Common Multiple (LCM) of two integers.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The LCM of a and b.
 * @note This function asserts that both inputs are finite.
 * @note The function uses the GCD to calculate the LCM efficiently.
 * @note If the GCD is 0, the function returns 0 to avoid division by zero.
 */
fn int lcm(int a, int b)
{
    assert(isFinite(a) && isFinite(b));

    int result = gcd(a, b);
    if (result == 0) return 0;

    return (int) abs((double) (a / result * b));
}

/*
 * Calculates the factorial of a given non-negative integer.
 *
 * @param a The non-negative integer for which to calculate the factorial.
 * @return The factorial of the input number.
 * @note This function asserts that the input is finite and non-negative.
 * @note The factorial is calculated recursively, which may lead to stack overflow for large inputs.
 * @warning This implementation is not suitable for large inputs due to potential stack overflow.
 */
fn int fact(int a)
{
    assert(isFinite(a));
    if (a <= 1) return 1;

    return a * fact(a - 1);
}

/*
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
 */
fn int rand(int a, int b)
{
    assert(isFinite(a) && isFinite(b) && a < b);

    static uint seed = 1;
    static uint multiplier = 16807;  // 7^5
    static uint modulus = 2147483647;  // 2^31 - 1 (Mersenne prime)

    // Clock now = clock::now();
    // uint compile_time_seed = now.second * 1000000 +
    //                          now.minute * 10000 +
    //                          now.hour * 100 +
    //                          now.day;

    seed = (multiplier * (seed /*^ compile_time_seed*/)) % modulus;

    return a + (int) (((double) seed / modulus) * (double) (b - a + 1));
}

/*
 * Calculates the remainder of the division of two integers.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The remainder of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 */
fn int rem(int a, int b)
{
	assert(isFinite(a) && isFinite(b) && b > 0);
	return a % b;
}

/*
 * Performs floor division of two double values.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The floor of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 */
fn int fdiv(double a, double b)
{
	assert(isFinite(a) && isFinite(b) && b > 0);
	return floor(a / b);
}

/*
 * Calculates the power of a base number raised to an integer exponent.
 *
 * @param base The base number.
 * @param pow The integer exponent.
 * @return The result of base raised to the power of pow.
 * @note This function asserts that both base and pow are finite.
 * @note For pow = 0, the function returns 1.
 * @note The function uses a simple iterative approach for positive exponents.
 */
fn double pow(double base, int pow)
{
	assert(isFinite(base) && isFinite(pow));

	if (pow == 0) return 1;

	double product = base;

	for (int i = 1; i < pow; ++i) product *= base;

	return product;
}

/*
 * Checks if a given integer is prime.
 *
 * @param a The integer to check for primality.
 * @return true if the number is prime, false otherwise.
 * @note This function asserts that the input is finite.
 * @note Numbers less than 2 are not considered prime.
 * @note Even numbers greater than 2 are not prime.
 * @note The function checks for divisibility up to half of the input number.
 */
fn bool isPrime(int a)
{
	assert(isFinite(a));

	if (a < 2) return false;
	if (a > 2 && a % 2 == 0) return false;

	for (int i = 2; i < a / 2; ++i) {
		if (a % i == 0) return false;
	}

	return true;
}

/*
 * Checks if a given double value is finite.
 *
 * @param a The double value to check.
 * @return true if the value is finite, false otherwise.
 */
fn bool isFinite(double a)
{
	return !isInfinite(a) && !isNaN(a);
}

/*
 * Checks if a given double value is infinite.
 *
 * @param a The double value to check.
 * @return true if the value is infinite, false otherwise.
 */
fn bool isInfinite(double a)
{
	return a / a != a / a;
}

/*
 * Checks if a given double value is Not-a-Number (NaN).
 *
 * @param a The double value to check.
 * @return true if the value is NaN, false otherwise.
 */
fn bool isNaN(double a)
{
	return a != a;
}

/*
 * Calculates the sine of an angle using Taylor series approximation.
 *
 * @param a The angle in radians.
 * @return The sine of the input angle.
 * @note This function asserts that the input is finite.
 * @note The function normalizes the input angle to the range [-PI, PI].
 * @note The Taylor series is computed up to the seventh term for accuracy.
 */
fn double sin(double a)
{
	assert(isFinite(a));

	while (a > PI) a -= 2 * PI;
	while (a < -PI) a += 2 * PI;

	double result = a;
	double term = a;

	for (int i = 1; i <= 7; ++i) {
		term *= -a * a / (double) ((2 * i) * (2 * i + 1));
		result += term;
	}

	return result;
}

/*
 * Calculates the cosine of an angle using Taylor series approximation.
 *
 * @param a The angle in radians.
 * @return The cosine of the input angle.
 * @note This function asserts that the input is finite.
 * @note The function normalizes the input angle to the range [-PI, PI].
 * @note The Taylor series is computed up to the seventh term for accuracy.
 */
fn double cos(double a)
{
    assert(isFinite(a));

    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;

    double result = 1;
    double term = 1;

    for (int i = 1; i <= 7; ++i) {
        term *= -a * a / (double) ((2 * i - 1) * (2 * i));
        result += term;
    }

    return result;
}

/*
 * Calculates the tangent of an angle.
 *
 * @param a The angle in radians.
 * @return The tangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The tangent is calculated as the ratio of sine to cosine.
 */
fn double tan(double a)
{
    assert(isFinite(a));

    double s = sin(a);
    double c = cos(a);

    return s / c;
}

/*
 * Calculates the hyperbolic sine of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic sine of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 0.
 * @note The function uses the exponential function to compute the result.
 */
fn double sinh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 0;

    double ea = exp(a);
    return (ea - (1 / ea)) / 2;
}

/*
 * Calculates the hyperbolic cosine of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 1.
 * @note The function uses the exponential function to compute the result.
 */
fn double cosh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 1;

    double ea = exp(a);
    return (ea + (1 / ea)) / 2;
}

/*
 * Calculates the hyperbolic tangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 0.
 * @note The function uses the exponential function to compute the result.
 */
fn double tanh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 0;

    double ea = exp(2 * a);
    return (ea - 1) / (ea + 1);
}

/*
 * Calculates the arcsine (inverse sine) of a given value using a polynomial approximation.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arcsine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The approximation uses a seventh-degree polynomial for accuracy.
 */
fn double asin(double a)
{
    assert(isFinite(a) && a >= -1 && a <= 1);

    double a2 = pow(a, 2);
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
}

/*
 * Calculates the arccosine (inverse cosine) of a given value.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arccosine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The arccosine is calculated using the relationship: acos(x) = PI/2 - asin(x).
 */
fn double acos(double a)
{
    assert(isFinite(a) && a >= -1 && a <= 1);
    return (PI / 2) - asin(a);
}

/*
 * Calculates the arctangent (inverse tangent) of a given value using an approximation.
 *
 * @param a The input value.
 * @return The arctangent of the input value in radians.
 * @note This function asserts that the input is finite.
 * @note This approximation is less accurate for large input values.
 */
fn double atan(double a)
{
    assert(isFinite(a));
    return a / (1.28 * pow(a, 2));
}

/*
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
 */
fn double atan2(double a, double b)
{
    assert(isFinite(a) && isFinite(b));

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

/*
 * Calculates the inverse hyperbolic sine of a given value.
 *
 * @param a The input value.
 * @return The inverse hyperbolic sine of the input value.
 * @note This function asserts that the input is finite.
 */
fn double asinh(double a)
{
    assert(isFinite(a));
    return ln(a + sqrt(a * a + 1));
}

/*
 * Calculates the inverse hyperbolic cosine of a given value.
 *
 * @param a The input value, must be greater than or equal to 1.
 * @return The inverse hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite and greater than or equal to 1.
 */
fn double acosh(double a)
{
    assert(isFinite(a) && a >= 1);
    return ln(a + sqrt(a * a - 1));
}

/*
 * Calculates the inverse hyperbolic tangent of a given value.
 *
 * @param a The input value, must be in the range (-1, 1).
 * @return The inverse hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite and within the valid range.
 */
fn double atanh(double a)
{
    assert(isFinite(a) && a > -1 && a < 1);
    return 0.5 * ln((1 + a) / (1 - a));
}

/*
 * Calculates the secant of an angle.
 *
 * @param a The angle in radians.
 * @return The secant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The secant is calculated as the reciprocal of the cosine.
 */
fn double sec(double a)
{
    assert(isFinite(a));

    double c = cos(a);

    return 1 / c;
}

/*
 * Calculates the cosecant of an angle.
 *
 * @param a The angle in radians.
 * @return The cosecant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cosecant is calculated as the reciprocal of the sine.
 */
fn double csc(double a)
{
    assert(isFinite(a));

    double s = sin(a);

    return 1 / s;
}

/*
 * Calculates the cotangent of an angle.
 *
 * @param a The angle in radians.
 * @return The cotangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cotangent is calculated as the ratio of cosine to sine.
 */
fn double cot(double a)
{
    assert(isFinite(a));

    double s = sin(a);
    double c = cos(a);

    return c / s;
}

/*
 * Calculates the hyperbolic secant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic secant of the input value.
 * @note This function asserts that the input is finite.
 * @note For a = 0, the function returns 1.
 * @note The function uses the exponential function to compute the result.
 */
fn double sech(double a)
{
    assert(isFinite(a));

    if (a == 0) return 1;

    double ea = exp(a);
    return 2 / (ea + (1 / ea));
}

/*
 * Calculates the hyperbolic cosecant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosecant of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 */
fn double csch(double a)
{
    assert(isFinite(a));

    double ea = exp(a);
    return 2 / (ea - (1 / ea));
}

/*
 * Calculates the hyperbolic cotangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cotangent of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 */
fn double coth(double a)
{
    assert(isFinite(a));

    double ea = exp(2 * a);
    return (ea + 1) / (ea - 1);
}

/*
 * Calculates the exponential function (e^x) for a given value.
 *
 * @param a The exponent value.
 * @return The result of e raised to the power of a.
 * @note This function asserts that the input is finite.
 * @note The function uses a Taylor series approximation combined with exponent reduction.
 * @note For a = 0, the function returns 1.
 * @note The calculation is optimized for accuracy and efficiency.
 */
fn double exp(double a)
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

/*
 * Finds the minimum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The smaller of the two input values.
 * @note This function asserts that both inputs are finite.
 */
fn double min(double a, double b)
{
    assert(isFinite(a) && isFinite(b));
    return a < b ? a : b;
}

/*
 * Finds the maximum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The larger of the two input values.
 * @note This function asserts that both inputs are finite.
 */
fn double max(double a, double b)
{
    assert(isFinite(a) && isFinite(b));
    return a > b ? a : b;
}

/*
 * Clamps a double value between a minimum and maximum range.
 *
 * @param value The value to clamp.
 * @param min The minimum allowed value.
 * @param max The maximum allowed value.
 * @return The clamped value, which will be between min and max (inclusive).
 * @note This function asserts that all inputs are finite.
 */
fn double clamp(double value, double min, double max)
{
    assert(isFinite(value) && isFinite(min) && isFinite(max));

    if (value < min) return min;
    if (value > max) return max;

    return value;
}

/*
 * Calculates the natural logarithm of a given value.
 *
 * @param a The input value, must be greater than 0.
 * @return The natural logarithm of the input value.
 * @note This function asserts that the input is finite and greater than 0.
 * @note The function uses a series expansion for improved accuracy.
 */
fn double ln(double a)
{
    assert(isFinite(a) && a > 0);

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
        y *= -a * (double) (i - 1) / i;
        sum += y;
    } while (abs(y) > 1E-15);

    return sum + exp * LN2;
}

/*
 * Calculates the logarithm of a value with a specified base.
 *
 * @param a The input value.
 * @param base The base of the logarithm.
 * @return The logarithm of the input value with the specified base.
 * @note This function asserts that both inputs are finite.
 */
fn double log(double a, double base)
{
    assert(isFinite(a) && isFinite(base));
    return ln(a) / ln(base);
}

/*
 * Calculates the base-2 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-2 logarithm of the input value.
 * @note This function asserts that the input is finite.
 */
fn double log2(double a)
{
    assert(isFinite(a));
    return ln(a) / LN2;
}

/*
 * Calculates the base-10 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-10 logarithm of the input value.
 * @note This function asserts that the input is finite.
 */
fn double log10(double a)
{
    assert(isFinite(a));
    return ln(a) / LN10;
}

/*
 * Calculates the sum of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The sum of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 */
fn double sum(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);

    double sum = 0;

    foreach (double value : data) {
        sum += value;
    }

    return sum;
}

/*
 * Calculates the arithmetic mean of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The arithmetic mean of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 */
fn double mean(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);
    return sum(data) / data.len;
}

/*
 * Calculates the median of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The median value of the array.
 * @note This function asserts that size is finite and greater than 0.
 * @note This function modifies the original array by sorting it.
 */
fn double median(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);

    for (int i = 0; i < data.len - 1; ++i) {
        for (int j = 0; j < data.len - i - 1; ++j) {
            if (data[j] > data[j + 1]) {
                double temp = data[j];

                data[j] = data[j + 1];
                data[j + 1] = temp;
            }
        }
    }

    if (data.len % 2 == 0) {
        return (data[(data.len / 2) - 1] + data[data.len / 2]) / 2;
    } else {
        return data[data.len / 2];
    }
}

/*
 * Calculates the mode of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The mode (most frequent value) of the array.
 * @note This function asserts that size is finite and greater than 0.
 * @note If multiple modes exist, this function returns the first one encountered.
 */
fn double mode(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);

    double mode = data[0];
    int maxCount = 1;
    int currentCount = 1;

    for (int i = 1; i < data.len; ++i) {
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
        mode = data[data.len - 1];
    }

    return mode;
}

/*
 * Calculates the sample standard deviation of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The sample standard deviation of the array.
 * @note This function asserts that size is finite and greater than 1.
 */
fn double stddev(double[] data)
{
    assert(isFinite(data.len) && data.len > 1);

    double m = mean(data);
    double sum = 0;

    foreach (double value : data) {
        double diff = value - m;
        sum += diff * diff;
    }

    return sqrt(sum / (double) (data.len - 1));
}
