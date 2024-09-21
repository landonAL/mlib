from typing import List
import struct

PI      = 3.1415926535897932
TAU     = 6.2831853071795864
E       = 2.7182818284590452
PHI     = 1.6180339887498948
LN2     = 0.6931471805599453
LN10    = 2.3025850929940457
LOG2E   = 1.4426950408889634
LOG10E  = 0.4342944819032518
EULER   = 0.5772156649015329
CATALAN = 0.9159655941772190

def to_bits(a):
    assert is_finite(a)

    s = struct.pack('>f', a)
    return struct.unpack('>l', s)[0]

def to_float(a):
    assert is_finite(a)

    s = struct.pack('>l', a)
    return struct.unpack('>f', s)[0]

# Converts degrees to radians.
#
# @param deg The angle in degrees to convert.
# @return The angle converted to radians.
# @note This function asserts that the input is finite.
def to_radian(deg: float) -> float:
    assert is_finite(deg)
    return deg * (PI / 180)

# Converts radians to degrees.
#
# @param rad The angle in radians to convert.
# @return The angle converted to degrees.
# @note This function asserts that the input is finite.
def to_degree(rad: float) -> float:
    assert is_finite(rad)
    return rad * (180 / PI)

# Calculates the floor of a given double value.
#
# @param a The double value to floor.
# @return The largest integer less than or equal to the input value.
# @note This function asserts that the input is finite.
def floor(a: float) -> int:
    assert is_finite(a)
    return int(a)

# Calculates the ceiling of a given double value.
#
# @param a The double value to calculate the ceiling for.
# @return The smallest integer greater than or equal to the input value.
# @note This function asserts that the input is finite.
def ceil(a: float) -> int:
    assert is_finite(a)
    return int(a + 1) if a > int(a) else int(a)

# Rounds a given double value to the nearest integer.
#
# @param a The double value to round.
# @return The nearest integer to the input value.
# @note This function asserts that the input is finite.
def round(a: float) -> int:
    assert is_finite(a)
    return int(a + 1) if a + 0.5 >= int(a + 1) else int(a)

# Calculates the absolute value of a given double.
#
# @param a The input double value.
# @return The absolute value of the input.
# @note This function asserts that the input is finite.
def abs(a: float) -> float:
    assert is_finite(a)
    return -a if a < 0 else a

# Calculates the square root of a given number using the Newton-Raphson method.
#
# @param a The number to calculate the square root of.
# @return The square root of the input number.
# @note This function asserts that the input is finite and non-negative.
def sqrt(a: float) -> float:
	assert is_finite(a) and a >= 0

	if a == 0: return 0

	b, root = a, 0

	for i in range(1, 18):
		b = root = (b + (a / b)) / 2

	return root

# Calculates the inverse square root of a given number.
#
# @param a The number to calculate the inverse square root of.
# @return The inverse square root of the input number.
# @note This function asserts that the input is finite.
def isqrt(a: float) -> float:
    assert is_finite(a)
    return 1 / sqrt(a)

# Calculates the quick inverse square root of a given number using the "Fast Inverse Square Root" algorithm.
#
# @param a The number to calculate the quick inverse square root of.
# @return An approximation of the inverse square root of the input number.
# @note This function uses a bit-level hack for initial guess and Newton's method for refinement.
# @note This function asserts that the input is finite.
def qisqrt(a: float) -> float:
    assert is_finite(a)

    x2 = a * 0.5
    y  = a
    i  = to_bits(y)
    i  = 0x5f3759df - (i >> 1)
    y  = to_float(i)
    y  *= (1.5 - (x2 * y * y))
    y  *= (1.5 - (x2 * y * y))

    return y

# Calculates the Greatest Common Divisor (GCD) of two integers using the Euclidean algorithm.
#
# @param a The first integer.
# @param b The second integer.
# @return The GCD of a and b.
# @note This function asserts that both inputs are finite.
# @note The function uses the absolute values of the inputs to handle negative numbers.
def gcd(a: int, b: int) -> int:
    assert is_finite(a) and is_finite(b)

    if a < 0: a = -a
    if b < 0: b = -b

    while b != 0:
        temp = b

        b = a % b
        a = temp

    return a

# Calculates the Least Common Multiple (LCM) of two integers.
#
# @param a The first integer.
# @param b The second integer.
# @return The LCM of a and b.
# @note This function asserts that both inputs are finite.
# @note The function uses the GCD to calculate the LCM efficiently.
# @note If the GCD is 0, the function returns 0 to avoid division by zero.
def lcm(a: int, b: int) -> int:
    assert is_finite(a) and is_finite(b)

    result = gcd(a, b)
    if result == 0: return 0

    return a // result * b if a * b >= 0 else -(a // result * b)

# Calculates the factorial of a given non-negative integer.
#
# @param a The non-negative integer for which to calculate the factorial.
# @return The factorial of the input number.
# @note This function asserts that the input is finite.
# @note The factorial is calculated iteratively to avoid stack overflow for large inputs.
def fact(a: int) -> int:
    assert is_finite(a)

    result = 1

    for i in range(2, a + 1):
        result *= i

    return result

# Generates a random integer within a specified range using a linear congruential generator (LCG).
#
# @param a The lower bound of the range (inclusive).
# @param b The upper bound of the range (inclusive).
# @return A random integer between a and b (inclusive).
# @note This function asserts that both inputs are finite and that a is less than b.
# @note The LCG uses a static seed, which is updated with each call to the function.
# @note The multiplier and modulus values are chosen to create a full-period generator.
def rand(a: int, b: int) -> int:
    assert is_finite(a) and is_finite(b) and a < b

    global seed
    seed = 1
    multiplier = 16807  # 7^5
    modulus = 2147483647  # 2^31 - 1 (Mersenne prime)

    seed = (multiplier * seed) % modulus

    return a + int(((seed / modulus) * (b - a + 1)))

# Calculates the remainder of the division of two integers.
#
# @param a The dividend.
# @param b The divisor.
# @return The remainder of a divided by b.
# @note This function asserts that both inputs are finite and that b is positive.
def rem(a: int, b: int) -> int:
	assert is_finite(a) and is_finite(b) and b > 0
	return a % b

# Performs floor division of two double values.
#
# @param a The dividend.
# @param b The divisor.
# @return The floor of a divided by b.
# @note This function asserts that both inputs are finite and that b is non-zero.
def fdiv(a: float, b: float) -> int:
	assert is_finite(a) and is_finite(b) and b > 0
	return floor(a / b)

# Calculates the power of a base number raised to an integer exponent.
#
# @param base The base number.
# @param pow The integer exponent.
# @return The result of base raised to the power of pow.
# @note This function asserts that both base and pow are finite.
# @note For pow = 0, the function returns 1.
# @note The function uses a simple iterative approach for positive exponents.
def pow(base: float, power: int) -> float:
    assert is_finite(base) and is_finite(power)

    if power == 0: return 1

    product = base

    for i in range(1, power): product *= base

    return product

# Checks if a given integer is prime.
#
# @param a The integer to check for primality.
# @return True if the number is prime, False otherwise.
# @note This function asserts that the input is finite.
# @note Numbers less than 2 are not considered prime.
# @note Even numbers greater than 2 are not prime.
# @note The function checks for divisibility up to half of the input number.
def is_prime(a: int) -> bool:
    assert is_finite(a)

    if a < 2: return False
    if a > 2 and a % 2 == 0: return False

    for i in range(2, a // 2):
        if a % i == 0: return False

    return True

# Checks if a given double value is finite.
#
# @param a The double value to check.
# @return True if the value is finite, False otherwise.
def is_finite(a: float) -> bool:
	return not is_infinite(a) and not is_nan(a)

# Checks if a given double value is infinite.
#
# @param a The double value to check.
# @return True if the value is infinite, False otherwise.
def is_infinite(a: float) -> bool:
	return a / a != a / a

# Checks if a given double value is Not-a-Number (NaN).
#
# @param a The double value to check.
# @return True if the value is NaN, False otherwise.
# @note This function compares the value to itself, as NaN is the only value that is not equal to itself.
def is_nan(a: float) -> bool:
	return a != a

# Calculates the sine of an angle using Taylor series approximation.
#
# @param a The angle in radians.
# @return The sine of the input angle.
# @note This function asserts that the input is finite.
# @note The function normalizes the input angle to the range [-PI, PI].
# @note The Taylor series is computed up to the 7th term for accuracy.
def sin(a: float) -> float:
    assert is_finite(a)

    while a > PI: a -= 2 * PI
    while a < -PI: a += 2 * PI

    result = term = a

    for i in range(1, 8):
        term *= -a * a / ((2 * i) * (2 * i + 1))
        result += term

    return result

# Calculates the cosine of an angle using Taylor series approximation.
#
# @param a The angle in radians.
# @return The cosine of the input angle.
# @note This function asserts that the input is finite.
# @note The function normalizes the input angle to the range [-PI, PI].
# @note The Taylor series is computed up to the 7th term for accuracy.
def cos(a: float) -> float:
    assert is_finite(a)

    while a > PI: a -= 2 * PI
    while a < -PI: a += 2 * PI

    result = term = 1

    for i in range(1, 8):
        term *= -a * a / ((2 * i - 1) * (2 * i))
        result += term

    return result

# Calculates the tangent of an angle.
#
# @param a The angle in radians.
# @return The tangent of the input angle.
# @note This function asserts that the input is finite.
# @note The tangent is calculated as the ratio of sine to cosine.
def tan(a: float) -> float:
    assert is_finite(a)

    s = sin(a)
    c = cos(a)

    return s / c

# Calculates the hyperbolic sine of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic sine of the input value.
# @note This function asserts that the input is finite.
# @note For a = 0, the function returns 0.
# @note The function uses the exponential function to compute the result.
def sinh(a: float) -> float:
    assert is_finite(a)

    if a == 0: return 0

    ea = exp(a)
    return (ea - (1 / ea)) / 2

# Calculates the hyperbolic cosine of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic cosine of the input value.
# @note This function asserts that the input is finite.
# @note For a = 0, the function returns 1.
# @note The function uses the exponential function to compute the result.
def cosh(a: float) -> float:
    assert is_finite(a)

    if a == 0: return 1

    ea = exp(a)
    return (ea + (1 / ea)) / 2

# Calculates the hyperbolic tangent of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic tangent of the input value.
# @note This function asserts that the input is finite.
# @note For a = 0, the function returns 0.
# @note The function uses the exponential function to compute the result.
def tanh(a: float) -> float:
    assert is_finite(a)

    if a == 0: return 0

    ea = exp(2 * a)
    return (ea - 1) / (ea + 1)

# Calculates the arcsine (inverse sine) of a given value using a polynomial approximation.
#
# @param a The input value, must be in the range [-1, 1].
# @return The arcsine of the input value in radians.
# @note This function asserts that the input is finite and within the valid range.
# @note The approximation uses a 7th-degree polynomial for accuracy.
def asin(a: float) -> float:
    assert is_finite(a) and a >= -1 and a <= 1

    a2 = pow(a, 2)
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)))

# Calculates the arccosine (inverse cosine) of a given value.
#
# @param a The input value, must be in the range [-1, 1].
# @return The arccosine of the input value in radians.
# @note This function asserts that the input is finite and within the valid range.
# @note The arccosine is calculated using the relationship: acos(x) = PI/2 - asin(x).
def acos(a: float) -> float:
    assert is_finite(a) and a >= -1 and a <= 1
    return (PI / 2) - asin(a)

# Calculates the arctangent (inverse tangent) of a given value using an approximation.
#
# @param a The input value.
# @return The arctangent of the input value in radians.
# @note This function asserts that the input is finite.
# @note This approximation is less accurate for large input values.
def atan(a: float) -> float:
    assert is_finite(a)
    return a / (1.28 * pow(a, 2))

# Calculates the arctangent of two variables (atan2).
#
# @param a The y-coordinate.
# @param b The x-coordinate.
# @return The angle in radians between the positive x-axis and the point (b, a).
# @note This function asserts that both inputs are finite.
# @note Special cases:
#       - If b is 0 and a > 0, returns PI/2
#       - If b is 0 and a < 0, returns -PI/2
#       - If b is 0 and a is 0, returns 0
#       - If b < 0, adjusts the result by adding or subtracting PI
def atan2(a: float, b: float) -> float:
    assert is_finite(a) and is_finite(b)

    if b == 0:
        if a > 0: return PI / 2
        if a < 0: return -PI / 2

        return 0

    result = atan(a / b)

    if b < 0:
        if a >= 0: return result + PI
        return result - PI

    return result

# Calculates the inverse hyperbolic sine of a given value.
#
# @param a The input value.
# @return The inverse hyperbolic sine of the input value.
# @note This function asserts that the input is finite.
def asinh(a: float) -> float:
    assert is_finite(a)
    return ln(a + sqrt(a * a + 1))

# Calculates the inverse hyperbolic cosine of a given value.
#
# @param a The input value, must be greater than or equal to 1.
# @return The inverse hyperbolic cosine of the input value.
# @note This function asserts that the input is finite and greater than or equal to 1.
def acosh(a: float) -> float:
    assert is_finite(a) and a >= 1
    return ln(a + sqrt(a * a - 1))

# Calculates the inverse hyperbolic tangent of a given value.
#
# @param a The input value, must be in the range (-1, 1).
# @return The inverse hyperbolic tangent of the input value.
# @note This function asserts that the input is finite and within the valid range.
def atanh(a: float) -> float:
    assert is_finite(a) and a > -1 and a < 1
    return 0.5 * ln((1 + a) / (1 - a))

# Calculates the secant of an angle.
#
# @param a The angle in radians.
# @return The secant of the input angle.
# @note This function asserts that the input is finite.
# @note The secant is calculated as the reciprocal of the cosine.
def sec(a: float) -> float:
    assert is_finite(a)

    c = cos(a)

    return 1 / c

# Calculates the cosecant of an angle.
#
# @param a The angle in radians.
# @return The cosecant of the input angle.
# @note This function asserts that the input is finite.
# @note The cosecant is calculated as the reciprocal of the sine.
def csc(a: float) -> float:
    assert is_finite(a)

    s = sin(a)

    return 1 / s

# Calculates the cotangent of an angle.
#
# @param a The angle in radians.
# @return The cotangent of the input angle.
# @note This function asserts that the input is finite.
# @note The cotangent is calculated as the ratio of cosine to sine.
def cot(a: float) -> float:
    assert is_finite(a)

    s = sin(a)
    c = cos(a)

    return c / s

# Calculates the hyperbolic secant of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic secant of the input value.
# @note This function asserts that the input is finite.
# @note For a = 0, the function returns 1.
# @note The function uses the exponential function to compute the result.
def sech(a: float) -> float:
    assert is_finite(a)

    if a == 0: return 1

    ea = exp(a)
    return 2 / (ea + (1 / ea))

# Calculates the hyperbolic cosecant of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic cosecant of the input value.
# @note This function asserts that the input is finite.
# @note The function uses the exponential function to compute the result.
def csch(a: float) -> float:
    assert is_finite(a)

    ea = exp(a)
    return 2 / (ea - (1 / ea))

# Calculates the hyperbolic cotangent of a given value.
#
# @param a The input value in radians.
# @return The hyperbolic cotangent of the input value.
# @note This function asserts that the input is finite.
# @note The function uses the exponential function to compute the result.
def coth(a: float) -> float:
    assert is_finite(a)

    ea = exp(2 * a)
    return (ea + 1) / (ea - 1)

# Calculates the exponential function (e^x) for a given value.
#
# @param a The exponent value.
# @return The result of e raised to the power of a.
# @note This function asserts that the input is finite.
# @note The function uses a Taylor series approximation combined with exponent reduction.
# @note For a = 0, the function returns 1.
# @note The calculation is optimized for accuracy and efficiency.
def exp(a: float) -> float:
    assert is_finite(a)

    if a == 0: return 1

    k = int(a * LOG2E)
    r = a - k * LN2
    result = r + 1
    term = r

    for i in range(2, 13):
        term *= r / i
        result += term

        if term < 1E-15 * result: break

    return result * pow(2, k)

# Finds the minimum of two double values.
#
# @param a The first double value to compare.
# @param b The second double value to compare.
# @return The smaller of the two input values.
# @note This function asserts that both inputs are finite.
def min(a: float, b: float) -> float:
    assert is_finite(a) and is_finite(b)
    return a if a < b else b

# Finds the maximum of two double values.
#
# @param a The first double value to compare.
# @param b The second double value to compare.
# @return The larger of the two input values.
# @note This function asserts that both inputs are finite.
def max(a: float, b: float) -> float:
    assert is_finite(a) and is_finite(b)
    return a if a > b else b

# Clamps a double value between a minimum and maximum range.
#
# @param value The value to clamp.
# @param min_val The minimum allowed value.
# @param max_val The maximum allowed value.
# @return The clamped value, which will be between min_val and max_val (inclusive).
# @note This function asserts that all inputs are finite.
def clamp(value: float, min_val: float, max_val: float) -> float:
    assert is_finite(value) and is_finite(min_val) and is_finite(max_val)

    if value < min_val: return min_val
    if value > max_val: return max_val

    return value

# Calculates the natural logarithm of a given value.
#
# @param a The input value, must be greater than 0.
# @return The natural logarithm of the input value.
# @note This function asserts that the input is finite and greater than 0.
# @note The function uses a series expansion for improved accuracy.
def ln(a: float) -> float:
    assert is_finite(a) and a > 0

    if a == 1: return 0

    exp = 0

    while a > 2:
        a /= 2
        exp += 1

    while a < 1:
        a *= 2
        exp -= 1

    a -= 1

    y = a
    sum = y
    i = 1

    while abs(y) > 1E-15:
        i += 1
        y *= -a * (i - 1) / i
        sum += y

    return sum + exp * LN2

# Calculates the logarithm of a value with a specified base.
#
# @param a The input value.
# @param base The base of the logarithm.
# @return The logarithm of the input value with the specified base.
# @note This function asserts that both inputs are finite.
def log(a: float, base: int) -> float:
    assert is_finite(a) and is_finite(base)
    return ln(a) / ln(base)

# Calculates the base-2 logarithm of a given value.
#
# @param a The input value.
# @return The base-2 logarithm of the input value.
# @note This function asserts that the input is finite.
def log2(a: float) -> float:
    assert is_finite(a)
    return ln(a) / LN2

# Calculates the base-10 logarithm of a given value.
#
# @param a The input value.
# @return The base-10 logarithm of the input value.
# @note This function asserts that the input is finite.
def log10(a: float) -> float:
    assert is_finite(a)
    return ln(a) / LN10

# Calculates the sum of an array of doubles.
#
# @param data List of double values.
# @return The sum of all elements in the array.
# @note This function asserts that the input list is finite and non-empty.
def sum(data: List[float]) -> float:
    assert is_finite(len(data)) and len(data) > 0

    sum = 0

    for i in range(len(data)):
        sum += data[i]

    return sum

# Calculates the arithmetic mean of an array of doubles.
#
# @param data List of double values.
# @return The arithmetic mean of all elements in the array.
# @note This function asserts that the input list is finite and non-empty.
def mean(data: List[float]) -> float:
    assert is_finite(len(data)) and len(data) > 0
    return sum(data) / len(data)

# Calculates the median of an array of doubles.
#
# @param data List of double values.
# @return The median value of the array.
# @note This function asserts that the input list is finite and non-empty.
# @note This function modifies the original list by sorting it.
def median(data: List[float]) -> float:
    assert is_finite(len(data)) and len(data) > 0

    size = len(data)

    for i in range(size - 1):
        for j in range(size - i - 1):
            if data[j] > data[j + 1]:
                data[j], data[j + 1] = data[j + 1], data[j]

    if size % 2 == 0: return (data[size // 2 - 1] + data[size // 2]) / 2
    else: return data[size // 2]

# Calculates the mode of an array of doubles.
#
# @param data List of double values.
# @return The mode (most frequent value) of the array.
# @note This function asserts that the input list is finite and non-empty.
# @note If multiple modes exist, this function returns the first one encountered.
def mode(data: List[float]) -> float:
    assert is_finite(len(data)) and len(data) > 0

    mode = data[0]
    max_count = 1
    current_count = 1

    for i in range(1, len(data)):
        if data[i] == data[i - 1]:
            current_count += 1
        else:
            if current_count > max_count:
                max_count = current_count
                mode = data[i - 1]

            current_count = 1

    if current_count > max_count:
        mode = data[-1]

    return mode

# Calculates the sample standard deviation of an array of doubles.
#
# @param data List of double values.
# @return The sample standard deviation of the array.
# @note This function asserts that the input list is finite and has more than one element.
def stddev(data: List[float]) -> float:
    assert is_finite(len(data)) and len(data) > 1

    size = len(data)
    m = mean(data)
    sum = 0

    for i in range(size):
        diff = data[i] - m
        sum += diff * diff

    return sqrt(sum / (size - 1))
