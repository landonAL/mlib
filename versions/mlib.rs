pub const PI:      f64 = 3.1415926535897932;
pub const TAU:     f64 = 6.2831853071795864;
pub const E:       f64 = 2.7182818284590452;
pub const PHI:     f64 = 1.6180339887498948;
pub const LN2:     f64 = 0.6931471805599453;
pub const LN10:    f64 = 2.3025850929940457;
pub const LOG2E:   f64 = 1.4426950408889634;
pub const LOG10E:  f64 = 0.4342944819032518;
pub const EULER:   f64 = 0.5772156649015329;
pub const CATALAN: f64 = 0.9159655941772190;

/**
 * Converts degrees to radians.
 *
 * @param deg The angle in degrees to convert.
 * @return The angle converted to radians.
 * @note This function asserts that the input is finite.
 **/
pub fn to_radian(deg: f64) -> f64 {
	assert!(is_finite(deg));
	return deg * (PI / 180.0);
}

/**
 * Converts radians to degrees.
 *
 * @param rad The angle in radians to convert.
 * @return The angle converted to degrees.
 * @note This function asserts that the input is finite.
 **/
pub fn to_degree(rad: f64) -> f64 {
	assert!(is_finite(rad));
	return rad * (180.0 / PI);
}

/**
 * Calculates the floor of a given double value.
 *
 * @param a The double value to floor.
 * @return The largest integer less than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
pub fn floor(a: f64) -> i32 {
    assert!(is_finite(a));
    return a as i32;
}

/**
 * Calculates the ceiling of a given double value.
 *
 * @param a The double value to calculate the ceiling for.
 * @return The smallest integer greater than or equal to the input value.
 * @note This function asserts that the input is finite.
 **/
pub fn ceil(a: f64) -> i32 {
	assert!(is_finite(a));
	return (a + 1) as i32;
}

/**
 * Rounds a given double value to the nearest integer.
 *
 * @param a The double value to round.
 * @return The nearest integer to the input value.
 * @note This function asserts that the input is finite.
 **/
pub fn round(a: f64) -> i32 {
	assert!(is_finite(a));
	return if a + 0.5 >= (a as i32 + 1) as f64 { a as i32 + 1 } else { a as i32 };
}

/**
 * Calculates the absolute value of a given double.
 *
 * @param a The input double value.
 * @return The absolute value of the input.
 * @note This function asserts that the input is finite.
 **/
pub fn abs(a: f64) -> f64 {
	assert!(is_finite(a));
	return if a < 0.0 { -a } else { a };
}

/**
 * Calculates the square root of a given number using the Newton-Raphson method.
 *
 * @param a The number to calculate the square root of.
 * @return The square root of the input number.
 * @note This function asserts that the input is finite and non-negative.
 **/
pub fn sqrt(a: f64) -> f64 {
	assert!(is_finite(a) && a >= 0.0);

	if a == 0.0 { return 0.0; }

	let mut b = a;
	let mut root = 0.0;

	for _ in 1..17 {
		b = (b + (a / b)) / 2.0;
		root = (b + (a / b)) / 2.0;
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
pub fn isqrt(a: f64) -> f64 {
    assert!(is_finite(a));
    return 1.0 / sqrt(a);
}

/**
 * Calculates the quick inverse square root of a given number using the "Fast Inverse Square Root" algorithm.
 *
 * @param a The number to calculate the quick inverse square root of.
 * @return An approximation of the inverse square root of the input number.
 * @note This function uses a bit-level hack for initial guess and Newton's method for refinement.
 * @note This function asserts that the input is finite.
 **/
pub fn qisqrt(a: f64) -> f64 {
    assert!(is_finite(a));

    let mut i: i64;
    let x2: f32;
    let mut y: f32;

    x2 = (a * 0.5) as f32;
    y = a as f32;
    i = unsafe { std::mem::transmute::<f32, i32>(y) as i64 };
    i = 0x5f3759df - (i >> 1);
    y = unsafe { std::mem::transmute::<i32, f32>(i as i32) };
    y *= 1.5 - (x2 * y * y);
    y *= 1.5 - (x2 * y * y);

    return y as f64;
}

/**
 * Calculates the nth root of a given number.
 *
 * @param a The number to calculate the root of.
 * @param exp The root exponent.
 * @return The nth root of the input number.
 * @note This function asserts that both inputs are finite, a is non-negative, and exp is positive.
 * @note For a = 0, the function returns 0.
 * @note The function uses the power operator for calculation.
 **/
pub fn root(a: f64, exp: i32) -> f64 {
    assert!(is_finite(a) && is_finite(exp.into()) && a >= 0.0 && exp > 0);
    if a == 0.0 { return 0.0; }

    let mut result = a;
    for _ in 0..16 {
        result = ((exp - 1) as f64 * result + a / pow(result, exp - 1)) / exp as f64;
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
pub fn gcd(a: i32, b: i32) -> i32 {
    assert!(is_finite(a.into()) && is_finite(b.into()));

    let mut a = abs(a as f64) as i32;
    let mut b = abs(b as f64) as i32;

    while b != 0 {
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
pub fn lcm(a: i32, b: i32) -> i32 {
    assert!(is_finite(a.into()) && is_finite(b.into()));

    let result = gcd(a, b);
    if result == 0 { return 0; }

    return abs((a / result * b) as f64) as i32;
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
 pub fn fact(a: i32) -> i32 {
     assert!(is_finite(a.into()));
     if a <= 1 { return 1; }

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
pub fn rand(a: i32, b: i32) -> i32 {
    assert!(is_finite(a.into()) && is_finite(b.into()) && a < b);

    static mut SEED: u32 = 1;
    const MULTIPLIER: u32 = 16807;  // 7^5
    const MODULUS: u32 = 2147483647;  // 2^31 - 1 (Mersenne prime)

    let compile_time_seed: u32 =
        (env!("CARGO_PKG_VERSION_MAJOR").parse::<u32>().unwrap_or(0) * 1000000) +
        (env!("CARGO_PKG_VERSION_MINOR").parse::<u32>().unwrap_or(0) * 100000) +
        (env!("CARGO_PKG_VERSION_PATCH").parse::<u32>().unwrap_or(0) * 10000) +
        (env!("CARGO_PKG_VERSION_PRE").chars().next().unwrap_or('0') as u32 - '0' as u32) * 1000 +
        (env!("CARGO_PKG_VERSION_PRE").chars().nth(1).unwrap_or('0') as u32 - '0' as u32) * 100 +
        (env!("CARGO_PKG_VERSION_PRE").chars().nth(2).unwrap_or('0') as u32 - '0' as u32) * 10;

    unsafe {
        SEED = (MULTIPLIER.wrapping_mul(SEED ^ compile_time_seed)) % MODULUS;
    }

    return a + ((unsafe { SEED } as f64 / MODULUS as f64) * (b - a + 1) as f64) as i32;
}

/**
 * Calculates the remainder of the division of two integers.
 *
 * @param a The dividend.
 * @param b The divisor.
 * @return The remainder of a divided by b.
 * @note This function asserts that both inputs are finite and that b is positive.
 **/
pub fn rem(a: i32, b: i32) -> i32 {
	assert!(is_finite(a.into()) && is_finite(b.into()) && b > 0);
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
pub fn fdiv(a: f64, b: f64) -> i32 {
	assert!(is_finite(a) && is_finite(b) && b > 0.0);
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
pub fn pow(base: f64, pow: i32) -> f64 {
	assert!(is_finite(base) && is_finite(pow.into()));

	if pow == 0 { return 1.0; }

	let mut product = base;

	for _ in 1..pow {
		product *= base;
	}

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
pub fn is_prime(a: i32) -> bool {
	assert!(is_finite(a.into()));

	if a < 2 { return false; }
	if a > 2 && a % 2 == 0 { return false;}

	for i in 2..a / 2 {
		if a % i == 0 { return false; }
	}

	return true;
}

/**
 * Checks if a given double value is finite.
 *
 * @param a The double value to check.
 * @return true if the value is finite, false otherwise.
 **/
pub fn is_finite(a: f64) -> bool {
	return !is_infinite(a) && !is_nan(a)
}

/**
 * Checks if a given double value is infinite.
 *
 * @param a The double value to check.
 * @return true if the value is infinite, false otherwise.
 **/
pub fn is_infinite(a: f64) -> bool {
	return a / a != a / a
}

/**
 * Checks if a given double value is Not-a-Number (NaN).
 *
 * @param a The double value to check.
 * @return true if the value is NaN, false otherwise.
 **/
pub fn is_nan(a: f64) -> bool {
	return a != a
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
pub fn sin(a: f64) -> f64 {
    assert!(is_finite(a));

    let mut a = a;

    while a > PI { a -= 2.0 * PI; }
    while a < -PI { a += 2.0 * PI; }

    let mut result = a;
    let mut term = a;

    for i in 1..=7 {
        term *= -a * a / ((2 * i) * (2 * i + 1)) as f64;
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
pub fn cos(a: f64) -> f64 {
    assert!(is_finite(a));

    let mut a = a;

    while a > PI { a -= 2.0 * PI; }
    while a < -PI { a += 2.0 * PI; }

    let mut result = 1.0;
    let mut term = 1.0;

    for i in 1..=7 {
        term *= -a * a / ((2 * i - 1) * (2 * i)) as f64;
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
pub fn tan(a: f64) -> f64 {
    assert!(is_finite(a));

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
pub fn sinh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 0.0; }

    let ea = exp(a);
    return (ea - (1.0 / ea)) / 2.0;
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
pub fn cosh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 1.0; }

    let ea = exp(a);
    return (ea + (1.0 / ea)) / 2.0;
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
pub fn tanh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 0.0; }

    let ea = exp(2.0 * a);
    return (ea - 1.0) / (ea + 1.0);
}

/**
 * Calculates the arcsine (inverse sine) of a given value using a polynomial approximation.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arcsine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The approximation uses a 7th-degree polynomial for accuracy.
 **/
pub fn asin(a: f64) -> f64 {
    assert!(is_finite(a) && a >= -1.0 && a <= 1.0);

    let a2 = pow(a, 2);
    return a + a * a2 * (1.0 / 6.0 + a2 * (3.0 / 40.0 + a2 * (5.0 / 112.0 + a2 * 35.0 / 1152.0)));
}

/**
 * Calculates the arccosine (inverse cosine) of a given value.
 *
 * @param a The input value, must be in the range [-1, 1].
 * @return The arccosine of the input value in radians.
 * @note This function asserts that the input is finite and within the valid range.
 * @note The arccosine is calculated using the relationship: acos(x) = PI/2 - asin(x).
 **/
pub fn acos(a: f64) -> f64 {
    assert!(is_finite(a) && a >= -1.0 && a <= 1.0);
    return (PI / 2.0) - asin(a);
}

/**
 * Calculates the arctangent (inverse tangent) of a given value using an approximation.
 *
 * @param a The input value.
 * @return The arctangent of the input value in radians.
 * @note This function asserts that the input is finite.
 * @note This approximation is less accurate for large input values.
 **/
pub fn atan(a: f64) -> f64 {
    assert!(is_finite(a));
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
pub fn atan2(a: f64, b: f64) -> f64 {
    assert!(is_finite(a) && is_finite(b));

    if b == 0.0 {
        if a > 0.0 { return PI / 2.0; }
        if a < 0.0 { return -PI / 2.0; }

        return 0.0;
    }

    let result = atan(a / b);

    if b < 0.0 {
        if a >= 0.0 { return result + PI; }
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
pub fn asinh(a: f64) -> f64 {
    assert!(is_finite(a));
    return ln(a + sqrt(a * a + 1.0));
}

/**
 * Calculates the inverse hyperbolic cosine of a given value.
 *
 * @param a The input value, must be greater than or equal to 1.
 * @return The inverse hyperbolic cosine of the input value.
 * @note This function asserts that the input is finite and greater than or equal to 1.
 **/
pub fn acosh(a: f64) -> f64 {
    assert!(is_finite(a) && a >= 1.0);
    return ln(a + sqrt(a * a - 1.0));
}

/**
 * Calculates the inverse hyperbolic tangent of a given value.
 *
 * @param a The input value, must be in the range (-1, 1).
 * @return The inverse hyperbolic tangent of the input value.
 * @note This function asserts that the input is finite and within the valid range.
 **/
pub fn atanh(a: f64) -> f64 {
    assert!(is_finite(a) && a > -1.0 && a < 1.0);
    return 0.5 * ln((1.0 + a) / (1.0 - a));
}

/**
 * Calculates the secant of an angle.
 *
 * @param a The angle in radians.
 * @return The secant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The secant is calculated as the reciprocal of the cosine.
 **/
pub fn sec(a: f64) -> f64 {
    assert!(is_finite(a));

    let c = cos(a);

    return 1.0 / cos(a);
}

/**
 * Calculates the cosecant of an angle.
 *
 * @param a The angle in radians.
 * @return The cosecant of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cosecant is calculated as the reciprocal of the sine.
 **/
pub fn csc(a: f64) -> f64 {
    assert!(is_finite(a));

    let s = sin(a);

    return 1.0 / sin(a);
}

/**
 * Calculates the cotangent of an angle.
 *
 * @param a The angle in radians.
 * @return The cotangent of the input angle.
 * @note This function asserts that the input is finite.
 * @note The cotangent is calculated as the ratio of cosine to sine.
 **/
pub fn cot(a: f64) -> f64 {
    assert!(is_finite(a));

    let s = sin(a);
    let c = cos(a);

    return cos(a) / sin(a);
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
pub fn sech(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 1.0; }

    let ea = exp(a);
    return 2.0 / (ea + (1.0 / ea));
}

/**
 * Calculates the hyperbolic cosecant of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cosecant of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
pub fn csch(a: f64) -> f64 {
    assert!(is_finite(a));

    let ea = exp(a);
    return 2.0 / (ea - (1.0 / ea));
}

/**
 * Calculates the hyperbolic cotangent of a given value.
 *
 * @param a The input value in radians.
 * @return The hyperbolic cotangent of the input value.
 * @note This function asserts that the input is finite.
 * @note The function uses the exponential function to compute the result.
 **/
pub fn coth(a: f64) -> f64 {
    assert!(is_finite(a));

    let ea = exp(2.0 * a);
    return (ea + 1.0) / (ea - 1.0);
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
pub fn exp(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 1.0; }

    let k = (a * LOG2E) as i32;
    let r = a - k as f64 * LN2;
    let mut result = r + 1.0;
    let mut term = r;

    for i in 2..=12 {
        term *= r / i as f64;
        result += term;

        if term < 1E-15 * result { break; }
    }

    return result * pow(2.0, k);
}

/**
 * Finds the minimum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The smaller of the two input values.
 * @note This function asserts that both inputs are finite.
 **/
pub fn min(a: f64, b: f64) -> f64 {
    assert!(is_finite(a) && is_finite(b));
    return if a < b { a } else { b };
}

/**
 * Finds the maximum of two double values.
 *
 * @param a The first double value to compare.
 * @param b The second double value to compare.
 * @return The larger of the two input values.
 * @note This function asserts that both inputs are finite.
 **/
pub fn max(a: f64, b: f64) -> f64 {
    assert!(is_finite(a) && is_finite(b));
    return if a > b { a } else { b };
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
pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
    assert!(is_finite(value) && is_finite(min) && is_finite(max));

    if value < min { return min; }
    if value > max { return max; }

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
pub fn ln(a: f64) -> f64 {
    assert!(is_finite(a) && a > 0.0);

    if a == 1.0 { return 0.0; }

    let mut exp = 0.0;
    let mut a = a;

    while a > 2.0 { a /= 2.0; exp += 1.0; }
    while a < 1.0 { a *= 2.0; exp -= 1.0; }

    a -= 1.0;

    let mut y = a;
    let mut sum = y;
    let mut i = 1;

    loop {
        i += 1;
        y *= -a * (i - 1) as f64 / i as f64;
        sum += y;

        if abs(y) <= 1E-15 { break; }
    }

    return sum + exp as f64 * LN2;
}

/**
 * Calculates the logarithm of a value with a specified base.
 *
 * @param a The input value.
 * @param base The base of the logarithm.
 * @return The logarithm of the input value with the specified base.
 * @note This function asserts that both inputs are finite.
 **/
pub fn log(a: f64, base: f64) -> f64 {
    assert!(is_finite(a) && is_finite(base));
    return ln(a) / ln(base);
}

/**
 * Calculates the base-2 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-2 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
pub fn log2(a: f64) -> f64 {
    assert!(is_finite(a));
    return ln(a) / LN2;
}

/**
 * Calculates the base-10 logarithm of a given value.
 *
 * @param a The input value.
 * @return The base-10 logarithm of the input value.
 * @note This function asserts that the input is finite.
 **/
pub fn log10(a: f64) -> f64 {
    assert!(is_finite(a));
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
pub fn sum(data: &[f64]) -> f64 {
    assert!(!data.is_empty());

    let mut sum = 0.0;

    for &value in data {
        sum += value;
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
pub fn mean(data: &[f64]) -> f64 {
    assert!(!data.is_empty());
    return sum(data) / data.len() as f64
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
pub fn median(data: &mut [f64]) -> f64 {
    assert!(!data.is_empty());

    data.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let size = data.len();
    if size % 2 == 0 { return (data[(size / 2) - 1] + data[size / 2]) / 2.0; } else { return data[size / 2]; }
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
pub fn mode(data: &[f64]) -> f64 {
    assert!(!data.is_empty());

    let mut sorted_data = data.to_vec();
    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut mode = sorted_data[0];
    let mut max_count = 1;
    let mut current_count = 1;

    for i in 1..sorted_data.len() {
        if sorted_data[i] == sorted_data[i - 1] {
            current_count += 1;
        } else {
            if current_count > max_count {
                max_count = current_count;
                mode = sorted_data[i - 1];
            }
            current_count = 1;
        }
    }

    if current_count > max_count { mode = sorted_data[sorted_data.len() - 1]; }

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
pub fn stddev(data: &[f64]) -> f64 {
    assert!(!data.is_empty() && data.len() > 1);

    let m = mean(data);
    let sum: f64 = data.iter().map(|&x| (x - m).powi(2)).sum();

    return sqrt(sum / (data.len() - 1) as f64);
}
