using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

// Note : When extending for full meshes (not only planes) this class would become about tetrahedral instead of triangles. 
public class FracTriangle {

    GameObject[] p;

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

    // Triangle property
    double halfArea;
    List<Vector<double>> forcesP = new List<Vector<double>>(); // Forces in the order of the points

    public FracTriangle(
        Matrix4x4 aWorld2MaterialTransform, 
        GameObject[] aParticleList, 
        double aMu,
        double aLambda,
        double aPhi,
        double aPsi)
    {
        world2Material = aWorld2MaterialTransform;

        p = aParticleList;
        Assert.AreEqual(3, aParticleList.Length, "FracTriangle: the particle list given have more or less than 3 points"); 

        mu = aMu;
        lambda = aLambda;
        phi = aPhi;
        psi = aPsi;
    }

    public bool Intersects(Vector3 planeNormal, Vector3 planePoint, float d, GameObject point, List<Vector3> intersectionTips)
    {
        bool intersects = false;

        GetSegmentPlaneIntersection(0, 1, planeNormal, d, intersectionTips);
        GetSegmentPlaneIntersection(1, 2, planeNormal, d, intersectionTips);
        GetSegmentPlaneIntersection(2, 0, planeNormal, d, intersectionTips);
        
        for(int i = 0; i < 3; ++i)
        {
            FracParticle fp = p[i].GetComponent<FracParticle>();

            if (!fp.SideIsSet)
            {
                // assign a side to points
                float side = Vector3.Dot(p[i].transform.position - planePoint, planeNormal);
                fp.Side = side > 0;
            }
            
        }

        return intersects; 
    }


    public void CalculateForces()
    {
        Update();

        Matrix<double> strain = Matrix<double>.Build.Dense(3, 3);
        Matrix<double> strainRate = Matrix<double>.Build.Dense(3, 3, 0);
        Matrix<double> elasticStress = Matrix<double>.Build.Dense(3, 3);
        Matrix<double> viscousStress = Matrix<double>.Build.Dense(3, 3);

        CalculateTensors(strain, strainRate, elasticStress, viscousStress);
        Matrix<double> stress = elasticStress + viscousStress;

        // For all points
        for (int i = 0; i < 3; ++i)
        {
            Vector<double> force = Vector<double>.Build.Dense(3);
            CalculatePointForce(force, i, stress);
            forcesP.Add(force);
        }

        // Separate stress tensor into tensile and compressive components
        CalculateTensilAndCompressive(stress);
    }

    private void CalculateTensors(Matrix<double> strainTensor, Matrix<double> strainRateTensor, 
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

    private void CalculateStrainTensor(Matrix<double> tensor, Vector<double> dxdu1, Vector<double> dxdu2, Vector<double> dxdu3)
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

    private void CalculateStrainRateTensor(Matrix<double> tensor, 
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

    private void CalculateStressTensor(Matrix<double> stress, Matrix<double> tensor, double a, double b)
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

    
    private void CalculatePointForce(Vector<double> outForce, int i, Matrix<double> totalStress)
    {
        outForce.Clear();

        Vector<double> pos = Vector<double>.Build.Dense(3);
        for (int comp = 0; comp < 3; ++comp)
            pos[comp] = p[i].transform.position[comp];

        // For all points
        for (int j = 0; j < 3; ++j)
        {
            double outerSum = 0;

            for (int k = 0; k < 3; ++k)
            {
                double innerSum = 0;

                for (int l = 0; l < 3; ++l)
                {
                    innerSum += beta[j, l] * beta[i, k] * totalStress[k, l];
                }

                outerSum += innerSum;
            }

            outForce += pos * outerSum;

        }
        
        outForce *= -halfArea;
    }

    private void CalculateTensilAndCompressive(Matrix<double> stress)
    {
        Evd<double> evd = stress.Evd();
        Matrix<double> tensileComp = Matrix<double>.Build.Dense(3, 3);
        Matrix<double> compressiveComp = Matrix<double>.Build.Dense(3, 3);

        for(int i = 0; i < 3; ++i)
        {
            double eval = evd.EigenValues[i].Real;
            Matrix<double> m = Utilities.BuildM(evd.EigenVectors.Column(i).Normalize(2)); 

            tensileComp += System.Math.Max(0, eval) * m;
            compressiveComp += System.Math.Min(0, eval) * m; 
        }

        for(int i = 0; i < 3; ++i)
        {
            Vector<double> ftensil = Vector<double>.Build.Dense(3);
            Vector<double> fcompressive = Vector<double>.Build.Dense(3);
            CalculatePointForce(ftensil, i, tensileComp);

            fcompressive = forcesP[i] - ftensil;

            FracParticle fp = p[i].GetComponent<FracParticle>();
            fp.AddTensilLoad(ftensil);
            fp.AddCompressiveLoad(fcompressive); 
        }

    }

    private void Update()
    {
        UpdateBeta();
        UpdateP();
        UpdateV();

        // Pre-compute area
        Vector3 ab = p[1].transform.position - p[0].transform.position;
        Vector3 ac = p[2].transform.position - p[0].transform.position;

        halfArea = Vector3.Cross(ab, ac).magnitude / 2.0d; 
        halfArea = System.Math.Abs(halfArea) / 2.0d;
    }

    private void UpdateBeta()
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

        Vector3 m1 = world2Material * p[0].transform.position;
        Vector3 m2 = world2Material * p[1].transform.position;
        Vector3 m3 = world2Material * p[2].transform.position;

        for (int i = 0; i < 3; ++i)
        {
            beta[i, 0] = m1[i];
            beta[i, 1] = m2[i];
            beta[i, 2] = m3[i];
        }

        beta = beta.Inverse();
    }

    private void UpdateP()
    {
        for (int i = 0; i < 3; ++i)
        {
            P[i, 0] = p[0].transform.position[i];
            P[i, 1] = p[1].transform.position[i];
            P[i, 2] = p[2].transform.position[i];
        }
    }

    private void UpdateV()
    {
        Rigidbody r1 = p[0].GetComponent<Rigidbody>();
        Rigidbody r2 = p[1].GetComponent<Rigidbody>();
        Rigidbody r3 = p[2].GetComponent<Rigidbody>();

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

    // This section uses the code from : http://stackoverflow.com/questions/3142469/determining-the-intersection-of-a-triangle-and-a-plane
    private float DistFromPlane(Vector3 p, Vector3 planeNormal, float d)
    {
        return Vector3.Dot(p, planeNormal) + d;
    }

    private void GetSegmentPlaneIntersection(int idA, int idB, Vector3 planeNormal, float d, List<Vector3> segTips)
    {
        Vector3 p0 = p[idA].transform.position;
        Vector3 p1 = p[idB].transform.position;

        float d1 = DistFromPlane(p0, planeNormal, d);
        float d2 = DistFromPlane(p1, planeNormal, d);

        float eps = 0.001f; 
        bool bP1OnPlane = Mathf.Abs(d1) < eps;
        bool bP2OnPlane = Mathf.Abs(d2) < eps;

        if (bP1OnPlane && !segTips.Exists(item => item.Equals(p0)))
            segTips.Add(p0);

        if (bP2OnPlane && !segTips.Exists(item => item.Equals(p1)))
            segTips.Add(p1);

        if (bP1OnPlane && bP2OnPlane)
            return;

        if (d1 * d2 > eps)
            return;

        float t = d1 / (d1 - d2);
        Vector3 res = p0 + t * (p1 - p0);
        if(!segTips.Exists(item => item.Equals(res)))
            segTips.Add(res); 
    }

}
