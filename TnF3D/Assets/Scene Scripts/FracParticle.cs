using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class FracParticle : MonoBehaviour {

    Vector<double> totalInternalForce = Vector<double>.Build.Dense(3, 0);

    public void ReinitializeForce()
    {
        totalInternalForce.Clear();
    }

    public void AddForce(Vector<double> aForce)
    {
        totalInternalForce += aForce;
    }
}
