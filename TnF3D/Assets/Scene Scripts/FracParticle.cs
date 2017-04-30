using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

public class FracParticle : MonoBehaviour {
    
    List<Vector<double>> tensilForces = new List<Vector<double>>();
    List<Vector<double>> compressiveForces = new List<Vector<double>>();

    Matrix<double> separationTensor = Matrix<double>.Build.Dense(3, 3);

    public void ReinitializeForce()
    {
        tensilForces.Clear();
        compressiveForces.Clear(); 
    }

    public void AddTensilLoad(Vector<double> aTensilForce)
    {
        tensilForces.Add(aTensilForce);
    }

    public void AddCompressiveLoad(Vector<double> aCompressiveForce)
    {
        compressiveForces.Add(aCompressiveForce);
    }

    public bool isFracturing(double toughness, List<Vector<double>> fracPlanesOrientation)
    {
        fracPlanesOrientation.Clear();

        CalculatePointSeparationTensor();
        Evd<double> evd = separationTensor.Evd();

        double vp = evd.EigenValues.AbsoluteMaximum().Real;
        bool fracturesAtNode = vp > toughness; // does at least one eigenvalue exceeds toughness
                                               //Debug.Log("vp = " + vp + "\n");

        if (vp > 0)
        {
            int i = 2 + 3;
            print(i); 
        }



        // Note : More efficient to simply loop on all since there's just 3 eigenvalue. Sorting is not that useful here.
        for (int i = 0; fracturesAtNode && i < evd.EigenValues.Count; ++i)
        {
            if (evd.EigenValues[i].Real > toughness)
            {
                Vector<double> fracPlaneOrientation = evd.EigenVectors.Column(i);
                fracPlanesOrientation.Add(fracPlaneOrientation);
            }
        }

        return fracturesAtNode;
    }

    private void CalculatePointSeparationTensor()
    {
        // There should be as many tensil components than compression components
        Assert.AreEqual(tensilForces.Count, compressiveForces.Count);

        Vector<double> f_plus = Vector<double>.Build.Dense(3, 0);
        Vector<double> f_min  = Vector<double>.Build.Dense(3, 0);
        
        Matrix<double> m_sumF_plus = Matrix<double>.Build.Dense(3, 3, 0);
        Matrix<double> m_sumF_min = Matrix<double>.Build.Dense(3, 3, 0);

        /* Compute sum(m(f)), f \in f+ or f \in f-
         * As well as calculate f+ & f- 
         */
        for (int i = 0; i < tensilForces.Count; ++i)
        {
            Vector<double> tf_i = tensilForces[i];
            Vector<double> cf_i = compressiveForces[i];

            m_sumF_plus += Utilities.BuildM(tf_i);
            m_sumF_min  += Utilities.BuildM(cf_i);

            f_plus += tf_i;
            f_min  += cf_i;
        }

        // Calculate m(f+) and m(f-)
        Matrix<double> mf_plus = Utilities.BuildM(f_plus);
        Matrix<double> mf_min = Utilities.BuildM(f_min);

        separationTensor = 1 / 2 * (-mf_plus + m_sumF_plus + mf_min - m_sumF_min);
    }

 
}
