using System;
using System.Diagnostics;

public static class mlib
{
    public const double PI      = 3.1415926535897932;
    public const double TAU     = 6.2831853071795864;
    public const double E       = 2.7182818284590452;
    public const double PHI     = 1.6180339887498948;
    public const double LN2     = 0.6931471805599453;
    public const double LN10    = 2.3025850929940457;
    public const double LOG2E   = 1.4426950408889634;
    public const double LOG10E  = 0.4342944819032518;
    public const double EULER   = 0.5772156649015329;
    public const double CATALAN = 0.9159655941772190;

    public static double ToRadian(double deg)
    {
    	Debug.Assert(IsFinite(deg));
    	return deg * (PI / 180);
    }

    public static double ToDegree(double rad)
    {
    	Debug.Assert(IsFinite(rad));
    	return rad * (180 / PI);
    }

    public static int Floor(double a)
    {
    	Debug.Assert(IsFinite(a));
    	return (int) a;
    }

    public static int Ceil(double a)
    {
    	Debug.Assert(IsFinite(a));
    	return a > (int) a ? (int) a + 1 : (int) a;
    }

    public static int Round(double a)
    {
        Debug.Assert(IsFinite(a));
        return a + 0.5 >= (int) a + 1
            ? (int) a + 1
            : (int) a;
    }

    public static double Abs(double a)
    {
        Debug.Assert(IsFinite(a));
        return a < 0 ? -a : a;
    }

    public static double Sqrt(double a)
    {
        Debug.Assert(IsFinite(a) && a >= 0);

        if (a == 0) return 0;

        double b = a, root = 0;

        for (int i = 1; i < 17; ++i)
        {
            b = root = (b + (a / b)) / 2;
        }

        return root;
    }

    public static double ISqrt(double a)
    {
        Debug.Assert(IsFinite(a));
        return 1 / Sqrt(a);
    }

    public static double QISqrt(double a)
    {
        Debug.Assert(IsFinite(a));

        long i;
        double x2, y;

        x2 = a * 0.5;
        y = a;
        i = BitConverter.DoubleToInt64Bits(y);
        i = 0x5f3759df - (i >> 1);
        y = BitConverter.Int64BitsToDouble(i);
        y *= (1.5 - (x2 * y * y));
        y *= (1.5 - (x2 * y * y));

        return y;
    }

    public static int GCD(int a, int b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b));

        if (a < 0) a = -a;
        if (b < 0) b = -b;

        while (b != 0)
        {
            int temp = b;

            b = a % b;
            a = temp;
        }

        return a;
    }

    public static int LCM(int a, int b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b));

        int result = GCD(a, b);
        if (result == 0) return 0;

        return a * b >= 0 ? a / result * b : -(a / result * b);
    }

    public static int Fact(int a)
    {
        Debug.Assert(IsFinite(a));

        int result = 1;

        for (int i = 2; i <= a; ++i)
        {
            result *= i;
        }

        return result;
    }

    public static int Rem(int a, int b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b) && b > 0);
        return a % b;
    }

    public static int Fdiv(double a, double b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b) && b > 0);
        return Floor(a / b);
    }

    public static double Pow(double baseV, int pow)
    {
        Debug.Assert(IsFinite(baseV) && IsFinite(pow));

        if (pow == 0) return 1;

        double product = baseV;

        for (int i = 1; i < pow; ++i) product *= baseV;

        return product;
    }

    public static bool IsPrime(int a)
    {
        Debug.Assert(IsFinite(a));

        if (a < 2) return false;
        if (a > 2 && a % 2 == 0) return false;

        for (int i = 2; i < a / 2; ++i)
        {
            if (a % i == 0) return false;
        }

        return true;
    }

    public static bool IsFinite(double a)
    {
        return !IsInfinite(a) && !IsNaN(a);
    }

    public static bool IsInfinite(double a)
    {
        return a / a != a / a;
    }

    public static bool IsNaN(double a)
    {
        return a != a;
    }

    public static double Sin(double a)
    {
        Debug.Assert(IsFinite(a));

        while (a > PI) a -= 2 * PI;
        while (a < -PI) a += 2 * PI;

        double result = a;
        double term = a;

        for (int i = 1; i <= 7; ++i)
        {
            term *= -a * a / ((2 * i) * (2 * i + 1));
            result += term;
        }

        return result;
    }

    public static double Cos(double a)
    {
        Debug.Assert(IsFinite(a));

        while (a > PI) a -= 2 * PI;
        while (a < -PI) a += 2 * PI;

        double result = 1;
        double term = 1;

        for (int i = 1; i <= 7; ++i)
        {
            term *= -a * a / ((2 * i - 1) * (2 * i));
            result += term;
        }

        return result;
    }

    public static double Tan(double a)
    {
        Debug.Assert(IsFinite(a));

        double s = Sin(a);
        double c = Cos(a);

        return s / c;
    }

    public static double Sinh(double a)
    {
        Debug.Assert(IsFinite(a));

        if (a == 0) return 0;

        double ea = Math.Exp(a);
        return (ea - (1 / ea)) / 2;
    }

    public static double Cosh(double a)
    {
        Debug.Assert(IsFinite(a));

        if (a == 0) return 1;

        double ea = Math.Exp(a);
        return (ea + (1 / ea)) / 2;
    }

    public static double Tanh(double a)
    {
        Debug.Assert(IsFinite(a));

        if (a == 0) return 0;

        double ea = Math.Exp(2 * a);
        return (ea - 1) / (ea + 1);
    }

    public static double Asin(double a)
    {
        Debug.Assert(IsFinite(a) && a >= -1 && a <= 1);

        double a2 = Pow(a, 2);
        return a + a * a2 * (1 / 6 + a2 * (3 / 40 + a2 * (5 / 112 + a2 * 35 / 1152)));
    }

    public static double Acos(double a)
    {
        Debug.Assert(IsFinite(a) && a >= -1 && a <= 1);
        return (PI / 2) - Asin(a);
    }

    public static double Atan(double a)
    {
        Debug.Assert(IsFinite(a));
        return a / (1.28 * Pow(a, 2));
    }

    public static double Atan2(double a, double b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b));

        if (b == 0)
        {
            if (a > 0) return PI / 2;
            if (a < 0) return -PI / 2;

            return 0;
        }

        double result = Atan(a / b);

        if (b < 0)
        {
            if (a >= 0) return result + PI;
            return result - PI;
        }

        return result;
    }

    public static double Asinh(double a)
    {
        Debug.Assert(IsFinite(a));
        return Math.Log(a + Math.Sqrt(a * a + 1));
    }

    public static double Acosh(double a)
    {
        Debug.Assert(IsFinite(a) && a >= 1);
        return Math.Log(a + Math.Sqrt(a * a - 1));
    }

    public static double Atanh(double a)
    {
        Debug.Assert(IsFinite(a) && a > -1 && a < 1);
        return 0.5 * Math.Log((1 + a) / (1 - a));
    }

    public static double Sec(double a)
    {
        Debug.Assert(IsFinite(a));

        double c = Cos(a);

        return 1 / c;
    }

    public static double Csc(double a)
    {
        Debug.Assert(IsFinite(a));

        double s = Sin(a);

        return 1 / s;
    }

    public static double Cot(double a)
    {
        Debug.Assert(IsFinite(a));

        double s = Sin(a);
        double c = Cos(a);

        return c / s;
    }

    public static double Sech(double a)
    {
        Debug.Assert(IsFinite(a));

        if (a == 0) return 1;

        double ea = Math.Exp(a);
        return 2 / (ea + (1 / ea));
    }

    public static double Csch(double a)
    {
        Debug.Assert(IsFinite(a));

        double ea = Math.Exp(a);
        return 2 / (ea - (1 / ea));
    }

    public static double Coth(double a)
    {
        Debug.Assert(IsFinite(a));

        double ea = Math.Exp(2 * a);
        return (ea + 1) / (ea - 1);
    }

    public static double Exp(double a)
    {
        Debug.Assert(IsFinite(a));

        if (a == 0) return 1;

        int k = (int) (a * LOG2E);
        double r = a - k * LN2;
        double result = r + 1;
        double term = r;

        for (int i = 2; i <= 12; ++i)
        {
            term *= r / i;
            result += term;

            if (term < 1E-15 * result) break;
        }

        return result * Pow(2, k);
    }

    public static double Min(double a, double b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b));
        return a < b ? a : b;
    }

    public static double Max(double a, double b)
    {
        Debug.Assert(IsFinite(a) && IsFinite(b));
        return a > b ? a : b;
    }

    public static double Clamp(double value, double min, double max)
    {
        Debug.Assert(IsFinite(value) && IsFinite(min) && IsFinite(max));

        if (value < min) return min;
        if (value > max) return max;

        return value;
    }

    public static double Ln(double a)
    {
        Debug.Assert(IsFinite(a) && a > 0);

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
        } while (Abs(y) > 1E-15);

        return sum + exp * LN2;
    }

    public static double Log(double a, double baseValue)
    {
        Debug.Assert(IsFinite(a) && IsFinite(baseValue));
        return Ln(a) / Ln(baseValue);
    }

    public static double Log2(double a)
    {
        Debug.Assert(IsFinite(a));
        return Ln(a) / LN2;
    }

    public static double Log10(double a)
    {
        Debug.Assert(IsFinite(a));
        return Ln(a) / LN10;
    }

    public static double Sum(double[] data)
    {
        Debug.Assert(data != null && data.Length > 0);

        double sum = 0;

        for (int i = 0; i < data.Length; ++i)
        {
            sum += data[i];
        }

        return sum;
    }

    public static double Mean(double[] data)
    {
        Debug.Assert(data != null && data.Length > 0);
        return Sum(data) / data.Length;
    }

    public static double Median(double[] data)
    {
        Debug.Assert(IsFinite(data.Length) && data.Length > 0);

        double[] sortedData = (double[]) data.Clone();
        Array.Sort(sortedData);

        int size = sortedData.Length;

        if (size % 2 == 0)
        {
            return (sortedData[(size / 2) - 1] + sortedData[size / 2]) / 2;
        }
        else
        {
            return sortedData[size / 2];
        }
    }

    public static double Mode(double[] data)
    {
        Debug.Assert(IsFinite(data.Length) && data.Length > 0);

        double mode = data[0];
        int maxCount = 1, currentCount = 1;

        for (int i = 1; i < data.Length; ++i)
        {
            if (data[i] == data[i - 1])
            {
                ++currentCount;
            }
            else
            {
                if (currentCount > maxCount)
                {
                    maxCount = currentCount;
                    mode = data[i - 1];
                }

                currentCount = 1;
            }
        }

        if (currentCount > maxCount)
        {
            mode = data[data.Length - 1];
        }

        return mode;
    }

    public static double StdDev(double[] data)
    {
        Debug.Assert(IsFinite(data.Length) && data.Length > 1);

        double m = Mean(data);
        double sum = 0;

        for (int i = 0; i < data.Length; i++)
        {
            double diff = data[i] - m;
            sum += diff * diff;
        }

        return Sqrt(sum / (data.Length - 1));
    }
}
