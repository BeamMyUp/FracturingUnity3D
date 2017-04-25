using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;
using MathNet.Numerics.LinearAlgebra;

// Note : When extending for full meshes (not only planes) this class would become about tetrahedral instead of triangles. 
public class FracTriangle {

    GameObject p1;
    GameObject p2;
    GameObject p3;

    Matrix4x4 world2Material; // The world to object space matrix of the mapped mesh

    // Physics data Matrices
    Matrix<double> beta = Matrix<double>.Build.Dense(3, 3, 0); // inverse material space positions matrix
    Matrix<double> P = Matrix <double>.Build.Dense(3, 3, 0); // world space positions matrix
    Matrix<double> V = Matrix<double>.Build.Dense(3, 3, 0); // world space velocities matrix

    // Material properties
    double mu; // Lame constant mu : material's rigidity
    double lambda; // Lame constant lambda : resistance to changes in volume (dilation)
    double phi;
    double psi; 

    public FracTriangle(
        Matrix4x4 aWorld2MaterialTransform, 
        GameObject aParticle1, 
        GameObject aParticle2, 
        GameObject aParticle3,
        double aMu,
        double aLambda,
        double aPhi,
        double aPsi)
    {
        world2Material = aWorld2MaterialTransform;

        p1 = aParticle1;
        p2 = aParticle2;
        p3 = aParticle3;

        mu = aMu;
        lambda = aLambda;
        phi = aPhi;
        psi = aPsi;

        UpdateMatrices();
        CalculateInternalStress(); 
    }

    public void CalculateInternalStress()
    {
        Matrix<double> strain = Matrix<double>.Build.Dense(3, 3);
        Matrix<double> strainRate = Matrix<double>.Build.Dense(3, 3, 0);
        Matrix<double> elasticStress = Matrix<double>.Build.Dense(3, 3);
        Matrix<double> viscousStress = Matrix<double>.Build.Dense(3, 3);

        CalculateTensors(strain, strainRate, elasticStress, viscousStress);
        Matrix<double> stress = elasticStress + viscousStress; 
    }

    public void CalculateTensors(Matrix<double> strainTensor, Matrix<double> strainRateTensor, 
                                 Matrix<double> elasticStressTensor, Matrix<double> viscousStressTensor)
    {
        // Calculate Partials
        Matrix<double> PB = P * beta;
        Matrix<double> VB = V * beta;

        /* dx/dui = P*Beta*Deltai, Deltai is 1 where i = j and 0 elsewhere
         * dx/dui = i-th column of P*Beta
         * 
         * dx'/dui = V*Beta*Deltai. For same reasons,
         * dx'/dui = i-th column of V*Beta
         */
        Vector<double> dxdu1 = PB.Column(0);
        Vector<double> dxdu2 = PB.Column(1);
        Vector<double> dxdu3 = PB.Column(2);

        Vector<double> dx_du1 = VB.Column(0);
        Vector<double> dx_du2 = VB.Column(1);
        Vector<double> dx_du3 = VB.Column(2);

        CalculateStrainTensor(strainTensor, dxdu1, dxdu2, dxdu3);
        CalculateStrainRateTensor(strainRateTensor, dxdu1, dxdu2, dxdu3, dx_du1, dx_du2, dx_du3);

        CalculateStressTensor(elasticStressTensor, strainTensor, lambda, mu);
        CalculateStressTensor(viscousStressTensor, strainRateTensor, phi, psi);
    }

    public void CalculateStrainTensor(Matrix<double> tensor, Vector<double> dxdu1, Vector<double> dxdu2, Vector<double> dxdu3)
    {
        // Diagonal
        tensor[0, 0] = dxdu1.DotProduct(dxdu1) - 1;
        tensor[1, 1] = dxdu2.DotProduct(dxdu2) - 1;
        tensor[2, 2] = dxdu3.DotProduct(dxdu3) - 1;

        // Symmetric sections
        tensor[0, 1] = tensor[1, 0] = dxdu1.DotProduct(dxdu2); 
        tensor[0, 2] = tensor[2, 0] = dxdu1.DotProduct(dxdu3); 
        tensor[1, 2] = tensor[2, 1] = dxdu2.DotProduct(dxdu3);
    }

    public void CalculateStrainRateTensor(Matrix<double> tensor, 
                                          Vector<double> dxdu1, Vector<double> dxdu2, Vector<double> dxdu3,
                                          Vector<double> dx_du1, Vector<double> dx_du2, Vector<double> dx_du3)
    {
        // Diagonal
        tensor[0, 0] = 2 * dxdu1.DotProduct(dx_du1);
        tensor[1, 1] = 2 * dxdu2.DotProduct(dx_du2);
        tensor[2, 2] = 2 * dxdu3.DotProduct(dx_du3);

        // Symmetric sections
        tensor[0, 1] = tensor[1, 0] = dxdu1.DotProduct(dx_du2) + dxdu2.DotProduct(dx_du1);
        tensor[0, 2] = tensor[2, 0] = dxdu1.DotProduct(dx_du3) + dxdu3.DotProduct(dx_du1);
        tensor[1, 2] = tensor[2, 1] = dxdu2.DotProduct(dx_du3) + dxdu3.DotProduct(dx_du2);
    }

    public void CalculateStressTensor(Matrix<double> stress, Matrix<double> tensor, double a, double b)
    {
        // TODO: optimise using the symmetry of the tensor (stress will also be symmetric)
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
            {
                int delta = i == j ? 1 : 0;
                double value = 0;

                for (int k = 0; k < 3; ++k)
                    value += (a * tensor[k, k] * delta) + (2 * b * tensor[i, j]);

                stress[i, j] = value;
            }
    }

    public void UpdateMatrices()
    {
        UpdateBeta();
        UpdateP();
        UpdateV();
    }

    public void UpdateBeta()
    {
        /* The Beta matrix is the inverse of the matrix made of the three points'
         * position in the material space
         * 
         * |m1.x  m2.x  m3.x|^-1
         * |m1.y  m2.y  m3.y|
         * |m1.z  m2.z  m3.z|
         * 
         * Where mi is a point in material space
         * Since all particles are GameObjects themselves, the material space is 
         * the Object Space of the mesh mapped onto the particles.
         */

        Vector3 m1 = world2Material * p1.transform.position;
        Vector3 m2 = world2Material * p2.transform.position;
        Vector3 m3 = world2Material * p3.transform.position;

        for (int i = 0; i < 3; ++i)
        {
            beta[i, 0] = m1[i];
            beta[i, 1] = m2[i];
            beta[i, 2] = m3[i];
        }

        beta = beta.Inverse();
    }

    public void UpdateP()
    {
        for (int i = 0; i < 3; ++i)
        {
            P[i, 0] = p1.transform.position[i];
            P[i, 1] = p2.transform.position[i];
            P[i, 2] = p3.transform.position[i];
        }
    }

    public void UpdateV()
    {
        Rigidbody r1 = p1.GetComponent<Rigidbody>();
        Rigidbody r2 = p2.GetComponent<Rigidbody>();
        Rigidbody r3 = p3.GetComponent<Rigidbody>();

        Assert.IsNotNull(r1, "FracTriangle.cs : P1 has no Rigidbody component");
        Assert.IsNotNull(r2, "FracTriangle.cs : P2 has no Rigidbody component");
        Assert.IsNotNull(r3, "FracTriangle.cs : P3 has no Rigidbody component");

        Vector3 v1 = r1.velocity;
        Vector3 v2 = r2.velocity;
        Vector3 v3 = r3.velocity;

        for (int i = 0; i < 3; ++i)
        {
            V[i, 0] = v1[i]; 
            V[i, 1] = v2[i]; 
            V[i, 2] = v3[i]; 
        }
    }

}
