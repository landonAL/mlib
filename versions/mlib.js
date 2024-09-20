const PI      = 3.1415926535897932;
const TAU     = 6.2831853071795864;
const E       = 2.7182818284590452;
const PHI     = 1.6180339887498948;
const LN2     = 0.6931471805599453;
const LN10    = 2.3025850929940457;
const LOG2E   = 1.4426950408889634;
const LOG10E  = 0.4342944819032518;
const EULER   = 0.5772156649015329;
const CATALAN = 0.9159655941772190;

function toRadian(deg) {
		if (!Number.isFinite(deg)) throw new Error("Input must be finite");
		return deg * (PI / 180);
}

function toDegree(rad) {
		if (!Number.isFinite(rad)) throw new Error("Input must be finite");
		return rad * (180 / PI);
}

function floor(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return Math.floor(a);
}

function ceil(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return a > Math.floor(a) ? Math.floor(a) + 1 : Math.floor(a);
}

function round(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return a + 0.5 >= Math.floor(a) + 1
				? Math.floor(a) + 1
				: Math.floor(a);
}

function abs(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return a < 0 ? -a : a;
}

function sqrt(a) {
		if (!Number.isFinite(a) || a < 0) throw new Error("Input must be finite and non-negative");

		if (a === 0) return 0;

		let b = a, root;

		for (let i = 1; i < 17; ++i) {
				b = root = (b + (a / b)) / 2;
		}

		return root;
}

function isqrt(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return 1 / sqrt(a);
}

function qisqrt(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let i;
    let x2 = a * 0.5;
    let y = a;

    let buffer = new ArrayBuffer(4);
    (new Float32Array(buffer))[0] = y;
    i = (new Uint32Array(buffer))[0];

    i = 0x5f3759df - (i >> 1);

    (new Uint32Array(buffer))[0] = i;
    y = (new Float32Array(buffer))[0];

    y *= (1.5 - (x2 * y * y));
    y *= (1.5 - (x2 * y * y));

    return y;
}

function gcd(a, b) {
    if (!Number.isFinite(a) || !Number.isFinite(b)) throw new Error("Inputs must be finite");

    a = abs(a);
    b = abs(b);

    while (b !== 0) {
        let temp = b;

        b = a % b;
        a = temp;
    }

    return a;
}

function lcm(a, b) {
    if (!Number.isFinite(a) || !Number.isFinite(b)) throw new Error("Inputs must be finite");

    let result = gcd(a, b);
    if (result === 0) return 0;

    return abs(a / result * b);
}

function fact(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");

		let result = 1;

		for (let i = 2; i <= a; ++i) {
				result *= i;
		}

		return result;
}

function rem(a, b) {
		if (!Number.isFinite(a) || !Number.isFinite(b) || b <= 0) throw new Error("Inputs must be finite and b must be positive");
		return a % b;
}

function fdiv(a, b) {
		if (!Number.isFinite(a) || !Number.isFinite(b) || b <= 0) throw new Error("Inputs must be finite and b must be positive");
		return floor(a / b);
}

function pow(base, pow) {
		if (!Number.isFinite(base) || !Number.isFinite(pow)) throw new Error("Inputs must be finite");

		if (pow === 0) return 1;

		let product = base;

		for (let i = 1; i < pow; ++i) product *= base;

		return product;
}

function isPrime(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");

		if (a < 2) return false;
		if (a > 2 && a % 2 === 0) return false;

		for (let i = 2; i < a / 2; ++i) {
				if (a % i === 0) return false;
		}

		return true;
}

function isFinite(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return !isInfinite(a) && !isNaN(a);
}

function isInfinite(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return a / a !== a / a;
}

function isNaN(a) {
		if (!Number.isFinite(a)) throw new Error("Input must be finite");
		return a !== a;
}

function sin(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;

    let result = a;
    let term = a;

    for (let i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i) * (2 * i + 1));
        result += term;
    }

    return result;
}

function cos(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;

    let result = 1;
    let term = 1;

    for (let i = 1; i <= 7; ++i) {
        term *= -a * a / ((2 * i - 1) * (2 * i));
        result += term;
    }

    return result;
}

function tan(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let s = sin(a);
    let c = cos(a);

    return s / c;
}

function sinh(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    if (a === 0) return 0;

    let ea = Math.exp(a);
    return (ea - (1 / ea)) / 2;
}

function cosh(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    if (a === 0) return 1;

    let ea = Math.exp(a);
    return (ea + (1 / ea)) / 2;
}

function tanh(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    if (a === 0) return 0;

    let ea = Math.exp(2 * a);
    return (ea - 1) / (ea + 1);
}

function asin(a) {
    if (!Number.isFinite(a) || a < -1 || a > 1) throw new Error("Input must be finite and between -1 and 1");

    let a2 = pow(a, 2);
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
}

function acos(a) {
    if (!Number.isFinite(a) || a < -1 || a > 1) throw new Error("Input must be finite and between -1 and 1");
    return (PI / 2) - asin(a);
}

function atan(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");
    return a / (1.28 * pow(a, 2));
}

function atan2(a, b) {
    if (!Number.isFinite(a) || !Number.isFinite(b)) throw new Error("Inputs must be finite");

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

function asinh(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");
    return Math.log(a + Math.sqrt(a * a + 1));
}

function acosh(a) {
    if (!Number.isFinite(a) || a < 1) throw new Error("Input must be finite and greater than or equal to 1");
    return Math.log(a + Math.sqrt(a * a - 1));
}

function atanh(a) {
    if (!Number.isFinite(a) || a <= -1 || a >= 1) throw new Error("Input must be finite and between -1 and 1 (exclusive)");
    return 0.5 * Math.log((1 + a) / (1 - a));
}

function sec(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let c = cos(a);

    return 1 / c;
}

function csc(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let s = sin(a);

    return 1 / s;
}

function cot(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let s = sin(a);
    let c = cos(a);

    return c / s;
}

function sech(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    if (a === 0) return 1;

    let ea = Math.exp(a);
    return 2 / (ea + (1 / ea));
}

function csch(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let ea = Math.exp(a);
    return 2 / (ea - (1 / ea));
}

function coth(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");

    let ea = Math.exp(2 * a);
    return (ea + 1) / (ea - 1);
}

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

        if (term < 1E-15 * result) break;
    }

    return result * pow(2, k);
}

function min(a, b) {
    if (!Number.isFinite(a) || !Number.isFinite(b)) throw new Error("Inputs must be finite");
    return a < b ? a : b;
}

function max(a, b) {
    if (!Number.isFinite(a) || !Number.isFinite(b)) throw new Error("Inputs must be finite");
    return a > b ? a : b;
}

function clamp(value, min, max) {
    if (!Number.isFinite(value) || !Number.isFinite(min) || !Number.isFinite(max)) throw new Error("Inputs must be finite");

    if (value < min) return min;
    if (value > max) return max;

    return value;
}

function ln(a) {
    if (!Number.isFinite(a) || a <= 0) throw new Error("Input must be finite and positive");

    if (a === 1) return 0;

    let exp = 0;

    while (a > 2) { a /= 2; exp++; }
    while (a < 1) { a *= 2; exp--; }

    a--;

    let y = a;
    let sum = y;
    let i = 1;

    do {
        i++;
        y *= -a * (i - 1) / i;
        sum += y;
    } while (abs(y) > 1E-15);

    return sum + exp * LN2;
}

function log(a, base) {
    if (!Number.isFinite(a) || !Number.isFinite(base)) throw new Error("Inputs must be finite");
    return ln(a) / ln(base);
}

function log2(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");
    return ln(a) / LN2;
}

function log10(a) {
    if (!Number.isFinite(a)) throw new Error("Input must be finite");
    return ln(a) / LN10;
}

function sum(data) {
    if (!Array.isArray(data) || data.length === 0) throw new Error("Input must be a non-empty array");

    let total = 0;

    for (let i = 0; i < data.length; ++i) {
        if (!Number.isFinite(data[i])) throw new Error("All elements must be finite numbers");
        total += data[i];
    }

    return total;
}

function mean(data) {
    if (!Array.isArray(data) || data.length === 0) throw new Error("Input must be a non-empty array");
    return sum(data) / data.length;
}
function median(data) {
    if (!Array.isArray(data) || data.length === 0) throw new Error("Input must be a non-empty array");

    for (let i = 0; i < data.length; i++) {
        if (!Number.isFinite(data[i])) throw new Error("All elements must be finite numbers");
    }

    let sortedData = [...data].sort((a, b) => a - b);
    let size = sortedData.length;

    if (size % 2 === 0) return (sortedData[(size / 2) - 1] + sortedData[size / 2]) / 2;
    else return sortedData[Math.floor(size / 2)];
}

function mode(data) {
    if (!Array.isArray(data) || data.length === 0) throw new Error("Input must be a non-empty array");

    for (let i = 0; i < data.length; i++) {
        if (!Number.isFinite(data[i])) throw new Error("All elements must be finite numbers");
    }

    let sortedData = [...data].sort((a, b) => a - b);
    let mode = sortedData[0];
    let maxCount = 1, currentCount = 1;

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

function stddev(data) {
    if (!Array.isArray(data) || data.length <= 1) throw new Error("Input must be an array with at least two elements");

    let m = mean(data);
    let sum = 0;

    for (let i = 0; i < data.length; i++) {
        if (!Number.isFinite(data[i])) throw new Error("All elements must be finite numbers");

        let diff = data[i] - m;
        sum += diff * diff;
    }

    return sqrt(sum / (data.length - 1));
}
