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

pub fn to_radian(deg: f64) -> f64 {
	assert!(is_finite(deg));
	return deg * (PI / 180.0);
}

pub fn to_degree(rad: f64) -> f64 {
	assert!(is_finite(rad));
	return rad * (180.0 / PI);
}

pub fn floor(a: f64) -> i32 {
    assert!(is_finite(a));
    return a.floor() as i32;
}

pub fn ceil(a: f64) -> i32 {
	assert!(is_finite(a));
	return a.ceil() as i32;
}

pub fn round(a: f64) -> i32 {
	assert!(is_finite(a));
	return if a + 0.5 >= (a as i32 + 1) as f64 { a as i32 + 1 } else { a as i32 };
}

pub fn abs(a: f64) -> f64 {
	assert!(is_finite(a));
	return if a < 0.0 { -a } else { a };
}

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

pub fn isqrt(a: f64) -> f64 {
    assert!(is_finite(a));
    return 1.0 / sqrt(a);
}

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

pub fn lcm(a: i32, b: i32) -> i32 {
    assert!(is_finite(a.into()) && is_finite(b.into()));

    let result = gcd(a, b);
    if result == 0 { return 0; }

    return abs((a / result * b) as f64) as i32;
}

pub fn fact(a: i32) -> i32 {
	assert!(is_finite(a.into()));

	let mut result = 1;

	for i in 2..=a {
		result *= i;
	}

	return result;
}

pub fn rem(a: i32, b: i32) -> i32 {
	assert!(is_finite(a.into()) && is_finite(b.into()) && b > 0);
	return a % b;
}

pub fn fdiv(a: f64, b: f64) -> i32 {
	assert!(is_finite(a) && is_finite(b) && b > 0.0);
	return floor(a / b);
}

pub fn pow(base: f64, pow: i32) -> f64 {
	assert!(is_finite(base) && is_finite(pow.into()));

	if pow == 0 { return 1.0; }

	let mut product = base;

	for _ in 1..pow {
		product *= base;
	}

	return product;
}

pub fn is_prime(a: i32) -> bool {
	assert!(is_finite(a.into()));

	if a < 2 { return false; }
	if a > 2 && a % 2 == 0 { return false;}

	for i in 2..a / 2 {
		if a % i == 0 { return false; }
	}

	return true;
}

pub fn is_finite(a: f64) -> bool {
	return !is_infinite(a) && !is_nan(a)
}

pub fn is_infinite(a: f64) -> bool {
	return a / a != a / a
}

pub fn is_nan(a: f64) -> bool {
	return a != a
}

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

pub fn tan(a: f64) -> f64 {
    assert!(is_finite(a));

    let s = sin(a);
    let c = cos(a);

    return s / c;
}

pub fn sinh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 0.0; }

    let ea = exp(a);
    return (ea - (1.0 / ea)) / 2.0;
}

pub fn cosh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 1.0; }

    let ea = exp(a);
    return (ea + (1.0 / ea)) / 2.0;
}

pub fn tanh(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 0.0; }

    let ea = exp(2.0 * a);
    return (ea - 1.0) / (ea + 1.0);
}

pub fn asin(a: f64) -> f64 {
    assert!(is_finite(a) && a >= -1.0 && a <= 1.0);

    let a2 = pow(a, 2);
    return a + a * a2 * (1.0 / 6.0 + a2 * (3.0 / 40.0 + a2 * (5.0 / 112.0 + a2 * 35.0 / 1152.0)));
}

pub fn acos(a: f64) -> f64 {
    assert!(is_finite(a) && a >= -1.0 && a <= 1.0);
    return (PI / 2.0) - asin(a);
}

pub fn atan(a: f64) -> f64 {
    assert!(is_finite(a));
    return a / (1.28 * pow(a, 2));
}

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

pub fn asinh(a: f64) -> f64 {
    assert!(is_finite(a));
    return ln(a + sqrt(a * a + 1.0));
}

pub fn acosh(a: f64) -> f64 {
    assert!(is_finite(a) && a >= 1.0);
    return ln(a + sqrt(a * a - 1.0));
}

pub fn atanh(a: f64) -> f64 {
    assert!(is_finite(a) && a > -1.0 && a < 1.0);
    return 0.5 * ln((1.0 + a) / (1.0 - a));
}

pub fn sec(a: f64) -> f64 {
    assert!(is_finite(a));

    let c = cos(a);

    return 1.0 / c;
}

pub fn csc(a: f64) -> f64 {
    assert!(is_finite(a));

    let s = sin(a);

    return 1.0 / s;
}

pub fn cot(a: f64) -> f64 {
    assert!(is_finite(a));

    let s = sin(a);
    let c = cos(a);

    return c / s;
}

pub fn sech(a: f64) -> f64 {
    assert!(is_finite(a));

    if a == 0.0 { return 1.0; }

    let ea = exp(a);
    return 2.0 / (ea + (1.0 / ea));
}

pub fn csch(a: f64) -> f64 {
    assert!(is_finite(a));

    let ea = exp(a);
    return 2.0 / (ea - (1.0 / ea));
}

pub fn coth(a: f64) -> f64 {
    assert!(is_finite(a));

    let ea = exp(2.0 * a);
    return (ea + 1.0) / (ea - 1.0);
}

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

pub fn min(a: f64, b: f64) -> f64 {
    assert!(is_finite(a) && is_finite(b));
    return if a < b { a } else { b };
}

pub fn max(a: f64, b: f64) -> f64 {
    assert!(is_finite(a) && is_finite(b));
    return if a > b { a } else { b };
}

pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
    assert!(is_finite(value) && is_finite(min) && is_finite(max));

    if value < min { return min; }
    if value > max { return max; }

    return value;
}

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

pub fn log(a: f64, base: f64) -> f64 {
    assert!(is_finite(a) && is_finite(base));
    return ln(a) / ln(base);
}

pub fn log2(a: f64) -> f64 {
    assert!(is_finite(a));
    return ln(a) / LN2;
}

pub fn log10(a: f64) -> f64 {
    assert!(is_finite(a));
    return ln(a) / LN10;
}

pub fn sum(data: &[f64]) -> f64 {
    assert!(!data.is_empty());

    let mut sum = 0.0;

    for &value in data {
        sum += value;
    }

    return sum;
}

pub fn mean(data: &[f64]) -> f64 {
    assert!(!data.is_empty());
    return sum(data) / data.len() as f64
}

pub fn median(data: &mut [f64]) -> f64 {
    assert!(!data.is_empty());

    data.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let size = data.len();
    if size % 2 == 0 { return (data[(size / 2) - 1] + data[size / 2]) / 2.0; } else { return data[size / 2]; }
}

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

pub fn stddev(data: &[f64]) -> f64 {
    assert!(!data.is_empty() && data.len() > 1);

    let m = mean(data);
    let sum: f64 = data.iter().map(|&x| (x - m).powi(2)).sum();

    return sqrt(sum / (data.len() - 1) as f64);
}
