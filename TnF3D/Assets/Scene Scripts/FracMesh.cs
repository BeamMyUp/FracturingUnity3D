using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;
using MathNet.Numerics.LinearAlgebra;

public class FracMesh : MonoBehaviour
{
    List<GameObject> particles = new List<GameObject>();
    List<FracTriangle> triangles = new List<FracTriangle>();

    // Material properties
    double mu; // Lame constant mu : material's rigidity
    double lambda; // Lame constant lambda : resistance to changes in volume (dilation)
    double phi;
    double psi;
    double tau;

    // Getters/Setters
    public int nParticles
    {
        get { return particles.Count; }
    }

    // TODO: find phi and psi meaning...
    public void InitializeFracMesh(double rigidity, double dilation, double aPhi, double aPsi, double toughness)
    {
        mu = rigidity;
        lambda = dilation;
        phi = aPhi;
        psi = aPsi;
        tau = toughness;
    }

    public void CreateMesh(int[] indices)
    {
        Vector3[] vertices = new Vector3[particles.Count];
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position;
        }

        // Create the mesh
        Mesh msh = new Mesh();
        msh.vertices = vertices;
        msh.triangles = indices;
        msh.RecalculateNormals();
        msh.RecalculateBounds();

        // Set up game object with mesh;
        gameObject.AddComponent(typeof(MeshRenderer));
        MeshFilter filter = gameObject.AddComponent(typeof(MeshFilter)) as MeshFilter;
        filter.mesh = msh;

        MeshRenderer mr = gameObject.GetComponent<MeshRenderer>();
        mr.material = new Material(Shader.Find("Diffuse"));
        mr.material.color = new Color(0.3f, 0.5f, 0.6f, 1.0f);

        // pin first particle
        particles[0].GetComponent<Rigidbody>().mass = 10;
    }

    public GameObject CreateParticle(Vector3 position)
    {
        GameObject particle = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        var fp = particle.AddComponent<FracParticle>();
        var rb = particle.AddComponent<Rigidbody>();
        rb.detectCollisions = true;
        rb.mass = 1;
        rb.drag = 5; 

        particle.transform.position = position;
        particle.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
        particles.Add(particle);

        fp.Initialize(particles.Count - 1);

        return particle;
    }

    public void CreateTriangle(int id1, int id2, int id3)
    {
        Matrix4x4 world2Object = gameObject.transform.worldToLocalMatrix;
        GameObject[] points = { particles[id1], particles[id2], particles[id3] };

        FracTriangle newTriangle = new FracTriangle(world2Object, points, mu, lambda, phi, psi);
        triangles.Add(newTriangle); 
    }

    public void AttachSpring(int iRb1, int iRb2)
    {
        // ---- Test on HingeJoint
        //HingeJoint joint = particles[iRb1].AddComponent<HingeJoint>();
        //joint.connectedBody = particles[iRb2].GetComponent<Rigidbody>();

        //JointSpring js = joint.spring;
        //js.spring = 2;
        //js.damper = 1;
        //js.targetPosition = 0;

        //joint.spring = js;
        //joint.useSpring = true;

        // ---- Test on Spring Joints
        //Rigidbody rb1 = particles[iRb1].GetComponent<Rigidbody>();
        //Rigidbody rb2 = particles[iRb2].GetComponent<Rigidbody>();
        //
        //SpringJoint joint = particles[iRb1].AddComponent<SpringJoint>();
        //joint.connectedBody = rb2;
        //
        //joint.anchor = rb1.position;
        //joint.connectedAnchor = joint.connectedBody.position;
        //
        //joint.spring = 20;
        //joint.damper = 2;
        //
        //joint.minDistance = 0;
        //joint.maxDistance = 0.0f;
        //joint.tolerance = 0; 

        // ---- Test on Custom Joints
        Spring joint = particles[iRb1].AddComponent<Spring>();
        joint.Initialize(particles[iRb1].GetComponent<Rigidbody>(), particles[iRb2].GetComponent<Rigidbody>());
    }

    private void UpdateMesh()
    {
        Vector3[] vertices = new Vector3[particles.Count];
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position;
        }

        Mesh m = gameObject.GetComponent<MeshFilter>().mesh;
        m.vertices = vertices;
        m.RecalculateNormals();
        m.RecalculateBounds();
    }

    private void Fracture()
    {
        // Reinitialize force accumulators and fracturing information of all particles
        for(int i = 0; i < particles.Count; ++i)
        {
            particles[i].GetComponent<FracParticle>().Reinitialize();
        }

        // Compute all forces applied to a particle and update force accumulators
        for(int i = 0; i < triangles.Count; ++i)
        {
            triangles[i].CalculateForces();
        }

        // Check if there are fractures at a particle's point
        for(int i = 0; i < particles.Count; ++i)
        {
            List<Vector<double>> planes = new List<Vector<double>>();
            FracParticle fpi = particles[i].GetComponent<FracParticle>();

            if (fpi.isFracturing(tau, planes))
            {
                Vector3 planeN = new Vector3();
                planeN.x = (float)planes[0].At(0);
                planeN.y = (float)planes[0].At(1);
                planeN.z = (float)planes[0].At(2);

                // Calculate distance from plane to origin
                Vector3 pos = particles[i].transform.position;
                float d = Mathf.Sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);

                for(int tri = 0; tri < triangles.Count; ++tri)
                {
                    // Check for intersection and assign side; only test on first plane for now
                    Vector3 tip = new Vector3();
                    GameObject particleOnPlane = null; 
                    bool isIntersectingElem = triangles[tri].Intersects(planeN, particles[i], d, tip, ref particleOnPlane);

                    if(isIntersectingElem)
                    {
                        CutElement(i, tip, triangles[tri], particleOnPlane);
                    }
                }
            }
        }

        // 5. Recalculate the mesh 
        Remeshing(); 
    }

    public void CutElement(int fracturePointId, Vector3 p1, FracTriangle triangle, GameObject particleOnPlane)
    {
        // Duplicate Fracture point
        GameObject originalP = particles[fracturePointId];
        FracParticle fpOriginal = originalP.GetComponent<FracParticle>();

        // Duplicate Particle and assign side of plane. 
        // The new particle is defined as always on the positive side of the plane
        GameObject duplicateP = DuplicateParticle(fracturePointId);
        FracParticle fpDupplicate = duplicateP.GetComponent<FracParticle>();

        fpOriginal.Side = false;
        fpDupplicate.Side = true;
        
        // If the other point is an already existing particle, then no other points to create
        GameObject secondPoint;
        if (particleOnPlane != null)
        {
            secondPoint = particleOnPlane;

            // HACK: This is a poor way of doing this search. Just trying to make it work here first.
            // Check if a spring between the originalPoint and the second point already exists
            bool alreadyExists = false; 
            for(int i = 0; i < fpOriginal.Springs.Count; ++i)
            {
                Spring s = fpOriginal.Springs[i];
                if ((s.P2.gameObject == secondPoint && s.P1.gameObject == originalP) ||
                    (s.P1.gameObject == secondPoint && s.P2.gameObject == originalP))
                {
                    alreadyExists = true;
                }
            }

            if(!alreadyExists)
            {
                Spring s = originalP.AddComponent<Spring>();
                s.Initialize(originalP.GetComponent<Rigidbody>(), secondPoint.GetComponent<Rigidbody>());
            }

            Spring ds = duplicateP.AddComponent<Spring>();
            ds.Initialize(duplicateP.GetComponent<Rigidbody>(), secondPoint.GetComponent<Rigidbody>());
        }

        else secondPoint = CreateParticle(p1);

        // Reassign springs
        for (int i = 0;  i < fpOriginal.Springs.Count; ++i)
        {
            Spring s = fpOriginal.Springs[i];
            GameObject otherPoint = s.GetOtherPoint(originalP);
            FracParticle sfp = otherPoint.GetComponent<FracParticle>();

            if(sfp.gameObject == secondPoint)
            {
                continue;
            }
            else
            {
                Assert.IsTrue(sfp.SideIsSet, "FracMesh : CutElement : Side is not set. Can't reassign springs");
            }

            // if spring particle doesn't have the same side as the initial particle, then the spring must be
            // re-assigned to the duplicate. Otherwise, nothing to change
            if (sfp.Side != fpOriginal.Side)
            {
                s.Initialize(otherPoint.GetComponent<Rigidbody>(), duplicateP.GetComponent<Rigidbody>());  
            }
        }
    }

    public void Remeshing()
    {
        List<int> indices = new List<int>();

        for (int i = 0; i < triangles.Count; ++i)
        {
            GameObject[] p = triangles[i].Points;
            FracParticle fp1 = p[0].GetComponent<FracParticle>();
            FracParticle fp2 = p[1].GetComponent<FracParticle>();
            FracParticle fp3 = p[2].GetComponent<FracParticle>();

            Vector3 A = p[0].transform.position;
            Vector3 B = p[1].transform.position;
            Vector3 C = p[2].transform.position;

            int p1 = fp1.Id;
            int p2 = fp2.Id;
            int p3 = fp3.Id;

            Vector3 n = new Vector3(0, 1, 0);
            float r = Vector3.Dot(n, Vector3.Cross(A - C, B - C));
            bool isClockwise = r > 0;

            if(!isClockwise)
            {
                p2 = fp3.Id;
                p3 = fp2.Id;
            }

            indices.Add(p1);
            indices.Add(p2);
            indices.Add(p3);
        }

        Vector3[] vertices = new Vector3[particles.Count];
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position;
        }

        // Create the mesh
        Mesh msh = new Mesh();
        msh.vertices = vertices;
        msh.triangles = indices.ToArray();
        msh.RecalculateNormals();
        msh.RecalculateBounds();
    }

    private GameObject DuplicateParticle(int idP2Copy)
    {
        GameObject iniP = particles[idP2Copy];
        Vector3 pos = iniP.transform.position;
        Vector3 vecPos = new Vector3(pos.x, pos.y, pos.z);

        return CreateParticle(vecPos);
    }

    void Update()
    {
        Fracture();
        UpdateMesh();
    }
}
