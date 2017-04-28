using MathNet.Numerics.LinearAlgebra;

class Utilities
{
    static public Matrix<double> BuildM(Vector<double> a)
    {
        Matrix<double> m = Matrix<double>.Build.Dense(a.Count, a.Count, 0);
        Matrix<double> vecA = Matrix<double>.Build.Dense(a.Count, 1);
        vecA.SetColumn(0, a);

        double norm = a.Norm(2);

        // if a != 0, m = aa^t/|a|
        if (norm != 0)
        {
            m = vecA * vecA.Transpose();
            m /= norm;
        }

        return m;
    }
}
