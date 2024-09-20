public class mlib {
    public static final double PI      = 3.1415926535897932;
    public static final double TAU     = 6.2831853071795864;
    public static final double E       = 2.7182818284590452;
    public static final double PHI     = 1.6180339887498948;
    public static final double LN2     = 0.6931471805599453;
    public static final double LN10    = 2.3025850929940457;
    public static final double LOG2E   = 1.4426950408889634;
    public static final double LOG10E  = 0.4342944819032518;
    public static final double EULER   = 0.5772156649015329;
    public static final double CATALAN = 0.9159655941772190;

    public static double toRadian(double deg) {
        assert isFinite(deg);
        return deg * (PI / 180);
    }

    public static double toDegree(double rad) {
        assert isFinite(rad);
        return rad * (180 / PI);
    }

    public static int floor(double a) {
        assert isFinite(a);
        return (int) a;
    }

    public static int ceil(double a) {
        assert isFinite(a);
        return a > (int) a ? (int) a + 1 : (int) a;
    }

    public static int round(double a) {
        assert isFinite(a);
        return a + 0.5 >= (int) a + 1
                ? (int) a + 1
                : (int) a;
    }

    public static double abs(double a) {
        assert isFinite(a);
        return a < 0 ? -a : a;
    }

    public static double sqrt(double a) {
        assert isFinite(a) && a >= 0;
        if (a == 0) return 0;

        double b = a, root;

        for (int i = 1; i < 17; ++i) {
            b = root = (b + (a / b)) / 2;
        }

        return root;
    }

    public static double isqrt(double a) {
        assert isFinite(a);
        return 1 / sqrt(a);
    }

    public static double qisqrt(double a) {
        assert isFinite(a);

        long i = 0;
        double x2 = 0, y = 0;

        x2 = a * 0.5;
        y = a;
        i = Double.floatToRawIntBits(y);
        i = 0x5f3759df - (i >> 1);
        y = Double.intBitsToFloat((int) i);
        y *= (1.5 - (x2 * y * y));
        y *= (1.5 - (x2 * y * y));

        return y;
    }

    public static int gcd(int a, int b) {
        assert isFinite(a) && isFinite(b);

        a = abs(a);
        b = abs(b);

        while (b != 0) {
            int temp = b;

            b = a % b;
            a = temp;
        }

        return a;
    }

    public static int lcm(int a, int b) {
        assert isFinite(a) && isFinite(b);

        int result = gcd(a, b);
        if (result == 0) return 0;

        return abs(a / result * b);
    }

    public static int fact(int a) {
        assert isFinite(a);

        int result = 1;

        for (int i = 2; i <= a; ++i) {
            result *= i;
        }

        return result;
    }

    public static int rem(int a, int b) {
        assert isFinite(a) && isFinite(b) && b > 0;
        return a % b;
    }

    public static int fdiv(double a, double b) {
        assert isFinite(a) && isFinite(b) && b > 0;
        return floor(a / b);
    }

    public static double pow(double base, int pow) {
        assert isFinite(base) && isFinite(pow);

        if (pow == 0) return 1;

        double product = base;

        for (int i = 1; i < pow; ++i) product *= base;

        return product;
    }

    public static boolean isPrime(int a) {
        assert isFinite(a);

        if (a < 2) return false;
        if (a > 2 && a % 2 == 0) return false;

        for (int i = 2; i < a / 2; ++i) {
            if (a % i == 0) return false;
        }

        return true;
    }

    public static boolean isFinite(double a) {
        return !isInfinite(a) && !isNaN(a);
    }

    public static boolean isInfinite(double a) {
        return a / a != a / a;
    }

    public static boolean isNaN(double a) {
        return a != a;
    }

    public static double sin(double a) {
        assert isFinite(a);

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

    public static double cos(double a) {
        assert isFinite(a);

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

    public static double tan(double a) {
        assert isFinite(a);

        double s = sin(a);
        double c = cos(a);

        return s / c;
    }

    public static double sinh(double a) {
        assert isFinite(a);

        if (a == 0) return 0;

        double ea = exp(a);
        return (ea - (1 / ea)) / 2;
    }

    public static double cosh(double a) {
        assert isFinite(a);

        if (a == 0) return 1;

        double ea = exp(a);
        return (ea + (1 / ea)) / 2;
    }

    public static double tanh(double a) {
        assert isFinite(a);

        if (a == 0) return 0;

        double ea = exp(2 * a);
        return (ea - 1) / (ea + 1);
    }

    public static double asin(double a) {
        assert isFinite(a) && a >= -1 && a <= 1;

        double a2 = pow(a, 2);
        return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
    }

    public static double acos(double a) {
        assert isFinite(a) && a >= -1 && a <= 1;
        return (PI / 2) - asin(a);
    }

    public static double atan(double a) {
        assert isFinite(a);
        return a / (1.28 * pow(a, 2));
    }

    public static double atan2(double a, double b) {
        assert isFinite(a) && isFinite(b);

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

    public static double asinh(double a) {
        assert isFinite(a);
        return ln(a + sqrt(a * a + 1));
    }

    public static double acosh(double a) {
        assert isFinite(a) && a >= 1;
        return ln(a + sqrt(a * a - 1));
    }

    public static double atanh(double a) {
        assert isFinite(a) && a > -1 && a < 1;
        return 0.5 * ln((1 + a) / (1 - a));
    }

    public static double sec(double a) {
        assert isFinite(a);

        double c = cos(a);

        return 1 / c;
    }

    public static double csc(double a) {
        assert isFinite(a);

        double s = sin(a);

        return 1 / s;
    }

    public static double cot(double a) {
        assert isFinite(a);

        double s = sin(a);
        double c = cos(a);

        return c / s;
    }

    public static double sech(double a) {
        assert isFinite(a);

        if (a == 0) return 1;

        double ea = exp(a);
        return 2 / (ea + (1 / ea));
    }

    public static double csch(double a) {
        assert isFinite(a);

        double ea = exp(a);
        return 2 / (ea - (1 / ea));
    }

    public static double coth(double a) {
        assert isFinite(a);

        if (a == 0) return isInfinite(a);

        double ea = exp(2 * a);
        return (ea + 1) / (ea - 1);
    }

    public static double exp(double a) {
        assert isFinite(a);

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

    public static double min(double a, double b) {
        assert isFinite(a) && isFinite(b);
        return a < b ? a : b;
    }

    public static double max(double a, double b) {
        assert isFinite(a) && isFinite(b);
        return a > b ? a : b;
    }

    public static double clamp(double value, double min, double max) {
        assert isFinite(value) && isFinite(min) && isFinite(max);

        if (value < min) return min;
        if (value > max) return max;

        return value;
    }

    public static double ln(double a) {
        assert isFinite(a) && a > 0;

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

    public static double log(double a, double base) {
        assert isFinite(a) && isFinite(base);
        return ln(a) / ln(base);
    }

    public static double log2(double a) {
        assert isFinite(a);
        return ln(a) / LN2;
    }

    public static double log10(double a) {
        assert isFinite(a);
        return ln(a) / LN10;
    }

    public static double sum(double[] data) {
        assert data != null && data.length > 0;

        double sum = 0;

        for (int i = 0; i < data.length; ++i) {
            sum += data[i];
        }

        return sum;
    }

    public static double mean(double[] data) {
        assert data != null && data.length > 0;
        return sum(data) / data.length;
    }

    public static double median(double[] data) {
        assert isFinite(data.length) && data.length > 0;

        for (int i = 0; i < data.length - 1; ++i) {
            for (int j = 0; j < data.length - i - 1; ++j) {
                if (data[j] > data[j + 1]) {
                    double temp = data[j];

                    data[j] = data[j + 1];
                    data[j + 1] = temp;
                }
            }
        }

        if (data.length % 2 == 0) {
            return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2;
        } else {
            return data[data.length / 2];
        }
    }

    public static double mode(double[] data) {
        assert data != null && data.length > 0;

        double mode = data[0];
        int maxCount = 1, currentCount = 1;

        for (int i = 1; i < data.length; ++i) {
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
            mode = data[data.length - 1];
        }

        return mode;
    }

    public static double stddev(double[] data) {
        assert data != null && data.length > 1;

        double m = mean(data);
        double sum = 0;

        for (int i = 0; i < data.length; i++) {
            double diff = data[i] - m;
            sum += diff * diff;
        }

        return sqrt(sum / (data.length - 1));
    }
}
