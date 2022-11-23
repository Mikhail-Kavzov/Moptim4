using System.Globalization;

internal class Program
{
    private static readonly int[] k_v = { 352000, 315000, 648000, 1027000, 242000 };
    private static readonly int[] sHalf = { 75, 130, 120, 100, 115 };
    private static double[] q = { 1, 1, 1, 1, 1 };

    private static void Main(string[] args)
    {
        const double e1 = 0.01;
        const double e2 = 0.02;
        const int m = 10;
        const int n = 5;
        bool isEnd = false;
        bool isPrev = false;
        int j = 0;
        int countIterations = 0;
        while (j < m && !isEnd)
        {
            int k = 0;
            bool isFirst = true;
            while (k <= n - 1)
            {
                ++countIterations;
                if (AbsGrad(q) < e1) //.7 |grad|<e1
                {
                    isEnd = true;
                    break;
                }
                var deriv = DfDq(q[k], sHalf[k], k_v[k]); //derivative
                var t0 = To(q[k], k_v[k], sHalf[k], deriv); // t0
                var prevF = F(q); //previous function
                var prevQ = q[k]; //previous q
                q[k] = RecalulateQ(q[k], deriv, t0); //.9
                if (Math.Abs(q[k] - prevQ) < e2 && Math.Abs(F(q) - prevF) < e2) //.10
                {
                    if (isPrev && isFirst) // j and j-1
                    {
                        isEnd = true;
                        break;
                    }
                    else
                    {
                        isPrev = true;
                    }
                }
                else
                {
                    isPrev = false;
                }
                isFirst = false;
                k++;
            }
            j++;
        }
        PrintResult(countIterations);
    }

    private static void PrintResult(int countIterations)
    {
        for (int i = 0; i < q.Length; i++)
        {
            Console.WriteLine($"q{i + 1} = {q[i]}");
        }
        Console.WriteLine($"F = {F(q)}");
        Console.WriteLine($"Count Iterations: {countIterations}");
    }

    private static double RecalulateQ(double q, double deriv, double t0) //9
    {
        return q - deriv * t0;
    }

    private static double To(double q, double kv, double sHalf, double deriv) //find t0
    {
        return (q - Math.Sqrt(kv / sHalf)) / deriv;
    }

    private static double AbsGrad(double[] q) //.7 |grad|
    {
        double result = 0;
        for (int i = 0; i < q.Length; i++)
        {
            var nq = DfDq(q[i], sHalf[i], k_v[i]);
            result += nq * nq;
        }
        return Math.Sqrt(result);
    }

    private static double DfDq(double q, double param1, double param2) //derivative
    {
        return param1 - param2 / (q * q);
    }

    private static double F(double[] q) //function (L)
    {
        double result = 0;
        for (int i = 0; i < q.Length; i++)
        {
            result += k_v[i] / q[i] + sHalf[i] * q[i];
        }
        return result;
    }
}