using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class FracParticle : MonoBehaviour {

    /* Let u = [u,v,w]^t be the position of the particle in the material's coordinate frame
     * In this case, this corresponds to local coordinates (object space)
     */

    Vector3 u;
    Matrix<double> strainTensor = Matrix<double>.Build.DenseIdentity(3, 3);
    Matrix<double> strainRateTensor = Matrix<double>.Build.Dense(3, 3, 0);
    Matrix4x4 prevL2W = new Matrix4x4(); 

    // Use this for initialization
    void Awake () {
        u = transform.localPosition;
        UpdateStrainTensor();
        UpdatePrevL2W(); 
	}

    public void UpdatePrevL2W()
    {
        Vector4 l2W0 = transform.localToWorldMatrix.GetColumn(0);
        Vector4 l2W1 = transform.localToWorldMatrix.GetColumn(1);
        Vector4 l2W2 = transform.localToWorldMatrix.GetColumn(2);
        Vector4 l2W3 = transform.localToWorldMatrix.GetColumn(3);

        prevL2W.SetColumn(0, l2W0);
        prevL2W.SetColumn(1, l2W1);
        prevL2W.SetColumn(2, l2W2);
        prevL2W.SetColumn(3, l2W3); 
    }

    public void UpdateStrainTensor()
    {
        // Pre-substract the Kronecker delta 
        strainTensor = Matrix<double>.Build.DenseIdentity(3, 3);
        strainTensor *= -1;

        Vector3 l2W0 = transform.localToWorldMatrix.GetColumn(0);
        Vector3 l2W1 = transform.localToWorldMatrix.GetColumn(1);
        Vector3 l2W2 = transform.localToWorldMatrix.GetColumn(2);

        // Only consider the 3x3 submatrix of the localToWorldMatrix 
        strainTensor[0, 0] += Vector3.Dot(l2W0, l2W0); 
        strainTensor[1, 1] += Vector3.Dot(l2W1, l2W1);
        strainTensor[2, 2] += Vector3.Dot(l2W2, l2W2);

        strainTensor[0, 1] = strainTensor[1, 0] = Vector3.Dot(l2W0, l2W1);
        strainTensor[0, 2] = strainTensor[2, 0] = Vector3.Dot(l2W0, l2W2);
        strainTensor[1, 2] = strainTensor[2, 1] = Vector3.Dot(l2W1, l2W2);
    }
	
    public void UpdateStrainRateTensor()
    {
        Rigidbody rb = gameObject.GetComponent<Rigidbody>();
        //Vector3 localOrigin = 
        //rb.GetRelativePointVelocity 
    }

    public void CalculateStress()
    {
        // calculate equation 7
        // calculate equation 8

        // add results
    }

    public void CalculatePotentialDensities()
    {
        // calculate elatic potential density (equation 9)
        // calculate damping potential density (equation 10)
    }

	// Update is called once per frame
	void Update () {
	    if(!prevL2W.Equals(transform.localToWorldMatrix))
        {
            //UpdatePrevL2W();
            //UpdateStrainTensor(); 
        }
	}
}
