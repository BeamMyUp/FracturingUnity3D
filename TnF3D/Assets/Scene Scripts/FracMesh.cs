using System.Collections;
using System.Collections.Generic;
using UnityEngine;
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
        particles[0].GetComponent<Rigidbody>().mass = 3;
    }

    public GameObject CreateParticle(Vector3 position)
    {
        GameObject particle = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        particle.AddComponent<FracParticle>();
        var rb = particle.AddComponent<Rigidbody>();
        rb.detectCollisions = true;
        rb.mass = 1;
        rb.drag = 5; 

        particle.transform.position = position;
        particle.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
        particles.Add(particle);

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
            if(particles[i].GetComponent<FracParticle>().isFracturing(tau, planes))
            {
                Vector3 planeN = new Vector3();
                planeN.x = (float)planes[0].At(0);
                planeN.y = (float)planes[0].At(1);
                planeN.z = (float)planes[0].At(2);

                Vector3 pos = particles[i].transform.position;
                float d = Mathf.Sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);

                // Duplicate Particle and add the new point at the end of the vector
                GameObject newP = DuplicateParticle(i);
                FracParticle fp = newP.GetComponent<FracParticle>();
                fp.Side = true; // the new particle is defined as on the positive side of the plane by convention

                for(int tri = 0; tri < triangles.Count; ++tri)
                {
                    // Check for intersection and assign side; only test on first plane for now
                    List<Vector3> tips = new List<Vector3>();
                    bool isIntersectingElem = triangles[tri].Intersects(planeN, pos, d, particles[i], tips);

                    if(isIntersectingElem)
                    {
                        // cut the element in two following tips point!! 
                    }
                }
            }
        }


        // 2. duplicate particles where the link breaks
        // 3. Re-evaluate links with the new particle
        // 4. apply force to the new particle
        // 5. Recalculate the mesh 
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
