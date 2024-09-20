module mlib;

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

fn double toRadian(double deg)
{
	assert(isFinite(deg));
	return deg * (PI / 180);
}

fn double toDegree(double rad)
{
    assert(isFinite(rad));
    return rad * (180 / PI);
}

fn int floor(double a)
{
	assert(isFinite(a));
	return (int) a;
}

fn int ceil(double a)
{
	assert(isFinite(a));
	return a > (int) a ? (int) a + 1 : (int) a;
}

fn int round(double a)
{
	assert(isFinite(a));
	return a + 0.5 >= (int) (a + 1)
					? (int) (a + 1)
					: (int) a;
}

fn double abs(double a)
{
	assert(isFinite(a));
	return a < 0 ? -a : a;
}

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

fn double isqrt(double a)
{
	assert(isFinite(a));
	return 1 / sqrt(a);
}

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

fn int lcm(int a, int b)
{
    assert(isFinite(a) && isFinite(b));

    int result = gcd(a, b);
    if (result == 0) return 0;

    return (int) abs((double) (a / result * b));
}

fn int fact(int a)
{
    assert(isFinite(a));

    int result = 1;

    for (int i = 2; i <= a; ++i)
    {
        result *= i;
    }

    return result;
}

fn int rem(int a, int b)
{
	assert(isFinite(a) && isFinite(b) && b > 0);
	return a % b;
}

fn int fdiv(double a, double b)
{
	assert(isFinite(a) && isFinite(b) && b > 0);
	return floor(a / b);
}

fn double pow(double base, int pow)
{
	assert(isFinite(base) && isFinite(pow));

	if (pow == 0) return 1;

	double product = base;

	for (int i = 1; i < pow; ++i) product *= base;

	return product;
}

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

fn bool isFinite(double a)
{
	return !isInfinite(a) && !isNaN(a);
}

fn bool isInfinite(double a)
{
	return a / a != a / a;
}

fn bool isNaN(double a)
{
	return a != a;
}

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

fn double tan(double a)
{
    assert(isFinite(a));

    double s = sin(a);
    double c = cos(a);

    return s / c;
}

fn double sinh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 0;

    double ea = exp(a);
    return (ea - (1 / ea)) / 2;
}

fn double cosh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 1;

    double ea = exp(a);
    return (ea + (1 / ea)) / 2;
}

fn double tanh(double a)
{
    assert(isFinite(a));

    if (a == 0) return 0;

    double ea = exp(2 * a);
    return (ea - 1) / (ea + 1);
}

fn double asin(double a)
{
    assert(isFinite(a) && a >= -1 && a <= 1);

    double a2 = pow(a, 2);
    return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
}

fn double acos(double a)
{
    assert(isFinite(a) && a >= -1 && a <= 1);
    return (PI / 2) - asin(a);
}

fn double atan(double a)
{
    assert(isFinite(a));
    return a / (1.28 * pow(a, 2));
}

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

fn double asinh(double a)
{
    assert(isFinite(a));
    return ln(a + sqrt(a * a + 1));
}

fn double acosh(double a)
{
    assert(isFinite(a) && a >= 1);
    return ln(a + sqrt(a * a - 1));
}

fn double atanh(double a)
{
    assert(isFinite(a) && a > -1 && a < 1);
    return 0.5 * ln((1 + a) / (1 - a));
}

fn double sec(double a)
{
    assert(isFinite(a));

    double c = cos(a);

    return 1 / c;
}

fn double csc(double a)
{
    assert(isFinite(a));

    double s = sin(a);

    return 1 / s;
}

fn double cot(double a)
{
    assert(isFinite(a));

    double s = sin(a);
    double c = cos(a);

    return c / s;
}

fn double sech(double a)
{
    assert(isFinite(a));

    if (a == 0) return 1;

    double ea = exp(a);
    return 2 / (ea + (1 / ea));
}

fn double csch(double a)
{
    assert(isFinite(a));

    double ea = exp(a);
    return 2 / (ea - (1 / ea));
}

fn double coth(double a)
{
    assert(isFinite(a));

    double ea = exp(2 * a);
    return (ea + 1) / (ea - 1);
}

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

fn double min(double a, double b)
{
    assert(isFinite(a) && isFinite(b));
    return a < b ? a : b;
}

fn double max(double a, double b)
{
    assert(isFinite(a) && isFinite(b));
    return a > b ? a : b;
}

fn double clamp(double value, double min, double max)
{
    assert(isFinite(value) && isFinite(min) && isFinite(max));

    if (value < min) return min;
    if (value > max) return max;

    return value;
}

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

fn double log(double a, double base)
{
    assert(isFinite(a) && isFinite(base));
    return ln(a) / ln(base);
}

fn double log2(double a)
{
    assert(isFinite(a));
    return ln(a) / LN2;
}

fn double log10(double a)
{
    assert(isFinite(a));
    return ln(a) / LN10;
}

fn double sum(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);

    double sum = 0;

    foreach (double value : data) {
        sum += value;
    }

    return sum;
}

fn double mean(double[] data)
{
    assert(isFinite(data.len) && data.len > 0);
    return sum(data) / data.len;
}

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
