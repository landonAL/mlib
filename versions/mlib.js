const PI = 3.1415926535897932;
const TAU = 6.2831853071795864;
const E = 2.7182818284590452;
const PHI = 1.6180339887498948;
const LN2 = 0.6931471805599453;
const LN10 = 2.3025850929940457;
const LOG2E = 1.4426950408889634;
const LOG10E = 0.4342944819032518;
const EULER = 0.5772156649015329;
const CATALAN = 0.915965594177219;

/**
 * Converts degrees to radians.
 *
 * @param deg The angle in degrees to convert.
 * @return The angle converted to radians.
 * @note This function asserts that the input is finite.
 **/
function toRadian(deg) {
  if (!Number.isFinite(deg)) throw new Error("Input must be finite");
  return deg * (PI / 180);
}

/**
 * Converts radians to degrees.
 *
 * @param rad The angle in radians to convert.
 * @return The angle converted to degrees.
 * @note This function asserts that the input is finite.
 **/
function toDegree(rad) {
  if (!Number.isFinite(rad)) throw new Error("Input must be finite");
  return rad * (180 / PI);
}

/**
 * Calculates the floor of a given double value.
 *
 * @param a The double value to floor.
 * @return The largest integer less than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
function floor(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  a = String.toString(a);
  return Number.parseInt(a);
}

/**
 * Calculates the ceiling of a given double value.
 *
 * @param a The double value to calculate the ceiling for.
 * @return The smallest integer greater than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
function ceil(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a > floor(a) ? floor(a) + 1 : floor(a);
}

/**
 * Rounds a given double value to the nearest integer.
 *
 * @param a The double value to round.
 * @return The nearest integer to the input value.
 * @note This function asserts that the input is finite.
 **/
function round(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a + 0.5 >= floor(a) + 1 ? floor(a) + 1 : floor(a);
}

/**
 * Calculates the absolute value of a given double.
 *
 * @param a The input double value.
 * @return The absolute value of the input.
 * @note This function asserts that the input is finite.
 **/
function abs(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a < 0 ? -a : a;
}

/**
 * Calculates the square root of a given number using the Newton-Raphson method.
 *
 * @param a The number to calculate the square root of.
 * @return The square root of the input number.
 * @note This function asserts that the input is finite and non-negative.
 **/
function sqrt(a) {
  if (!Number.isFinite(a) || a < 0)
    throw new Error("Input must be finite and non-negative");

  if (a === 0) return 0;

  let b = a,
    root;

  for (let i = 1; i < 17; ++i) {
    b = root = (b + a / b) / 2;
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
function isqrt(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return 1 / sqrt(a);
}

/**
 * Calculates the quick inverse square root of a given number using the "Fast Inverse Square Root" algorithm.
 *
 * @param a The number to calculate the quick inverse square root of.
 * @return An approximation of the inverse square root of the input number.
 * @note This function uses a bit-level hack for initial guess and Newton's method for refinement.
 * @note This function asserts that the input is finite.
 **/
function qisqrt(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let i;
  let x2 = a * 0.5;
  let y = a;

  let buffer = new ArrayBuffer(4);
  new Float32Array(buffer)[0] = y;
  i = new Uint32Array(buffer)[0];

  i = 0x5f3759df - (i >> 1);

  new Uint32Array(buffer)[0] = i;
  y = new Float32Array(buffer)[0];

  y *= 1.5 - x2 * y * y;
  y *= 1.5 - x2 * y * y;

  return y;
}

/**
 * Calculates the nth root of a given number.
 *
 * @param {number} a The number to calculate the root of.
 * @param {number} exp The root exponent.
 * @return {number} The nth root of the input number.
 * @throws {Error} If inputs are not finite, a is negative, or exp is not positive.
 * @note For a = 0, the function returns 0.
 * @note The function uses the power operator for calculation.
 */
function root(a, exp) {
  if (!Number.isFinite(a) || !Number.isFinite(exp) || a < 0 || exp <= 0) {
    throw new Error(
      "Invalid input: a must be non-negative, exp must be positive, and both must be finite",
    );
  }

  if (a === 0) return 0;

  let result = a;
  for (let i = 0; i < 16; i++) {
    result = ((exp - 1) * result + a / pow(result, exp - 1)) / exp;
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
function gcd(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b))
    throw new Error("Inputs must be finite");

  a = abs(a);
  b = abs(b);

  while (b !== 0) {
    let temp = b;

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
function lcm(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b))
    throw new Error("Inputs must be finite");

  let result = gcd(a, b);
  if (result === 0) return 0;

  return abs((a / result) * b);
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
function fact(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  if (a <= 1) return 1;

  return a * fact(a - 1);
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
function rand(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b) || a >= b)
    throw new Error("Inputs must be finite and a must be less than b");

  let seed = 1;
  const multiplier = 16807; // 7^5
  const modulus = 2147483647; // 2^31 - 1 (Mersenne prime)

  const compile_time_seed = new Date().getTime() % modulus;

  seed = (multiplier * (seed ^ compile_time_seed)) % modulus;

  return a + Math.floor((seed / modulus) * (b - a + 1));
}

/**
 * Calculates the remainder of the division of two integers.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The remainder of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 **/
function rem(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b) || b <= 0)
    throw new Error("Inputs must be finite and b must be positive");
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
function fdiv(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b) || b <= 0)
    throw new Error("Inputs must be finite and b must be positive");
  return floor(a / b);
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
function pow(base, pow) {
  if (!Number.isFinite(base) || !Number.isFinite(pow))
    throw new Error("Inputs must be finite");

  if (pow === 0) return 1;

  let product = base;

  for (let i = 1; i < pow; ++i) product *= base;

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
function isPrime(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a < 2) return false;
  if (a > 2 && a % 2 === 0) return false;

  for (let i = 2; i < a / 2; ++i) {
    if (a % i === 0) return false;
  }

  return true;
}

/**
 * Checks if a given double value is finite.
 *
 * @param a The double value to check.
 * @return true if the value is finite, false otherwise.
 **/
function isFinite(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return !isInfinite(a) && !isNaN(a);
}

/**
 * Checks if a given double value is infinite.
 *
 * @param a The double value to check.
 * @return true if the value is infinite, false otherwise.
 **/
function isInfinite(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a / a !== a / a;
}

/**
 * Checks if a given double value is Not-a-Number (NaN).
 *
 * @param a The double value to check.
 * @return true if the value is NaN, false otherwise.
 **/
function isNaN(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a !== a;
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
function sin(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  while (a > PI) a -= 2 * PI;
  while (a < -PI) a += 2 * PI;

  let result = a;
  let term = a;

  for (let i = 1; i <= 7; ++i) {
    term *= (-a * a) / (2 * i * (2 * i + 1));
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
function cos(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  while (a > PI) a -= 2 * PI;
  while (a < -PI) a += 2 * PI;

  let result = 1;
  let term = 1;

  for (let i = 1; i <= 7; ++i) {
    term *= (-a * a) / ((2 * i - 1) * (2 * i));
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
function tan(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let s = sin(a);
  let c = cos(a);

  return s / c;
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
function sinh(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a === 0) return 0;

  let ea = exp(a);
  return (ea - 1 / ea) / 2;
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
function cosh(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a === 0) return 1;

  let ea = exp(a);
  return (ea + 1 / ea) / 2;
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
function tanh(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a === 0) return 0;

  let ea = exp(2 * a);
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
function asin(a) {
  if (!Number.isFinite(a) || a < -1 || a > 1)
    throw new Error("Input must be finite and between -1 and 1");

  let a2 = pow(a, 2);
  return (
    a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + (a2 * 35) / 1152)))
  );
}

/**
 * Calculates the arccosine (inverse cosine) of a given value.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arccosine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The arccosine is calculated using the relationship: acos(x) = PI/2 - asin(x).
 **/
function acos(a) {
  if (!Number.isFinite(a) || a < -1 || a > 1)
    throw new Error("Input must be finite and between -1 and 1");
  return PI / 2 - asin(a);
}

/**
 * Calculates the arctangent (inverse tangent) of a given value using an approximation.
 *
 * @param a The input value.
 * @return The arctangent of the input value in radians.
 * @note This function asserts that the input is finite.
 * @note This approximation is less accurate for large input values.
 **/
function atan(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return a / (1.28 * pow(a, 2));
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
function atan2(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b))
    throw new Error("Inputs must be finite");

  if (b === 0) {
    if (a > 0) return PI / 2;
    if (a < 0) return -PI / 2;

    return 0;
  }

  let result = atan(a / b);

  if (b < 0) {
    if (a >= 0) return result + PI;
    return result - PI;
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
function asinh(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return log(a + sqrt(a * a + 1));
}

/**
 * Calculates the inverse hyperbolic cosine of a given value.
 *
 * @param a The input value, must be greater than or equal to 1.
 * @return The inverse hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite and greater than or equal to 1.
 **/
function acosh(a) {
  if (!Number.isFinite(a) || a < 1)
    throw new Error("Input must be finite and greater than or equal to 1");
  return log(a + sqrt(a * a - 1));
}

/**
 * Calculates the inverse hyperbolic tangent of a given value.
 *
 * @param a The input value, must be in the range (-1, 1).
 * @return The inverse hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite and within the valid range.
 **/
function atanh(a) {
  if (!Number.isFinite(a) || a <= -1 || a >= 1)
    throw new Error("Input must be finite and between -1 and 1 (exclusive)");
  return 0.5 * log((1 + a) / (1 - a));
}

/**
 * Calculates the secant of an angle.
 *
 * @param a The angle in radians.
 * @return The secant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The secant is calculated as the reciprocal of the cosine.
 **/
function sec(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let c = cos(a);

  return 1 / c;
}

/**
 * Calculates the cosecant of an angle.
 *
 * @param a The angle in radians.
 * @return The cosecant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cosecant is calculated as the reciprocal of the sine.
 **/
function csc(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let s = sin(a);

  return 1 / s;
}

/**
 * Calculates the cotangent of an angle.
 *
 * @param a The angle in radians.
 * @return The cotangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cotangent is calculated as the ratio of cosine to sine.
 **/
function cot(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let s = sin(a);
  let c = cos(a);

  return c / s;
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
function sech(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a === 0) return 1;

  let ea = exp(a);
  return 2 / (ea + 1 / ea);
}

/**
 * Calculates the hyperbolic cosecant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosecant of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
function csch(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let ea = exp(a);
  return 2 / (ea - 1 / ea);
}

/**
 * Calculates the hyperbolic cotangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cotangent of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
function coth(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  let ea = exp(2 * a);
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
function exp(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");

  if (a === 0) return 1;

  let k = floor(a * LOG2E);
  let r = a - k * LN2;
  let result = r + 1;
  let term = r;

  for (let i = 2; i <= 12; ++i) {
    term *= r / i;
    result += term;

    if (term < 1e-15 * result) break;
  }

  return result * pow(2, k);
}

/**
 * Finds the minimum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The smaller of the two input values.
 * @note This function asserts that both inputs are finite.
 **/
function min(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b))
    throw new Error("Inputs must be finite");
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
function max(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b))
    throw new Error("Inputs must be finite");
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
function clamp(value, min, max) {
  if (!Number.isFinite(value) || !Number.isFinite(min) || !Number.isFinite(max))
    throw new Error("Inputs must be finite");

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
function ln(a) {
  if (!Number.isFinite(a) || a <= 0)
    throw new Error("Input must be finite and positive");

  if (a === 1) return 0;

  let exp = 0;

  while (a > 2) {
    a /= 2;
    exp++;
  }
  while (a < 1) {
    a *= 2;
    exp--;
  }

  a--;

  let y = a;
  let sum = y;
  let i = 1;

  do {
    i++;
    y *= (-a * (i - 1)) / i;
    sum += y;
  } while (abs(y) > 1e-15);

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
function log(a, base) {
  if (!Number.isFinite(a) || !Number.isFinite(base))
    throw new Error("Inputs must be finite");
  return ln(a) / ln(base);
}

/**
 * Calculates the base-2 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-2 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
function log2(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return ln(a) / LN2;
}

/**
 * Calculates the base-10 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-10 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
function log10(a) {
  if (!Number.isFinite(a)) throw new Error("Input must be finite");
  return ln(a) / LN10;
}

/**
 * Calculates the sum of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The sum of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 **/
function sum(data) {
  if (!Array.isArray(data) || data.length === 0)
    throw new Error("Input must be a non-empty array");

  let total = 0;

  for (let i = 0; i < data.length; ++i) {
    if (!Number.isFinite(data[i]))
      throw new Error("All elements must be finite numbers");
    total += data[i];
  }

  return total;
}

/**
 * Calculates the arithmetic mean of an array of doubles.
 *
 * @param data Pointer to the array of doubles.
 * @param size The number of elements in the array.
 * @return The arithmetic mean of all elements in the array.
 * @note This function asserts that size is finite and greater than 0.
 **/
function mean(data) {
  if (!Array.isArray(data) || data.length === 0)
    throw new Error("Input must be a non-empty array");
  return sum(data) / data.length;
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
function median(data) {
  if (!Array.isArray(data) || data.length === 0)
    throw new Error("Input must be a non-empty array");

  for (let i = 0; i < data.length; i++) {
    if (!Number.isFinite(data[i]))
      throw new Error("All elements must be finite numbers");
  }

  let sortedData = [...data].sort((a, b) => a - b);
  let size = sortedData.length;

  if (size % 2 === 0)
    return (sortedData[size / 2 - 1] + sortedData[size / 2]) / 2;
  else return sortedData[floor(size / 2)];
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
function mode(data) {
  if (!Array.isArray(data) || data.length === 0)
    throw new Error("Input must be a non-empty array");

  for (let i = 0; i < data.length; i++) {
    if (!Number.isFinite(data[i]))
      throw new Error("All elements must be finite numbers");
  }

  let sortedData = [...data].sort((a, b) => a - b);
  let mode = sortedData[0];
  let maxCount = 1,
    currentCount = 1;

  for (let i = 1; i < sortedData.length; ++i) {
    if (sortedData[i] === sortedData[i - 1]) {
      ++currentCount;
    } else {
      if (currentCount > maxCount) {
        maxCount = currentCount;
        mode = sortedData[i - 1];
      }

      currentCount = 1;
    }
  }

  if (currentCount > maxCount) mode = sortedData[sortedData.length - 1];

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
function stddev(data) {
  if (!Array.isArray(data) || data.length <= 1)
    throw new Error("Input must be an array with at least two elements");

  let m = mean(data);
  let sum = 0;

  for (let i = 0; i < data.length; i++) {
    if (!Number.isFinite(data[i]))
      throw new Error("All elements must be finite numbers");

    let diff = data[i] - m;
    sum += diff * diff;
  }

  return sqrt(sum / (data.length - 1));
}
