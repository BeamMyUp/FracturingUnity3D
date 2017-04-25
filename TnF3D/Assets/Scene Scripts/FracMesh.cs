using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FracMesh
{
    GameObject go = new GameObject("Fracturable Mesh");
    List<GameObject> particles = new List<GameObject>();
    List<FracTriangle> triangles = new List<FracTriangle>();

    // Material properties
    double mu; // Lame constant mu : material's rigidity
    double lambda; // Lame constant lambda : resistance to changes in volume (dilation)
    double phi;
    double psi;

    // TODO: find phi and psi meaning...
    public FracMesh(double rigidity, double dilation, double aPhi, double aPsi)
    {
        mu = rigidity;
        lambda = dilation;
        phi = aPhi;
        psi = aPsi;
    }

    public int nParticles
    {
        get { return particles.Count; }
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
        go.AddComponent(typeof(MeshRenderer));
        MeshFilter filter = go.AddComponent(typeof(MeshFilter)) as MeshFilter;
        filter.mesh = msh;

        MeshRenderer mr = go.GetComponent<MeshRenderer>();
        mr.material = new Material(Shader.Find("Diffuse"));
        mr.material.color = new Color(0.3f, 0.5f, 0.6f, 1.0f);
    }

    public void UpdateMesh()
    {
        Vector3[] vertices = new Vector3[particles.Count]; 
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position;
        }

        Mesh m = go.GetComponent<MeshFilter>().mesh;
        m.vertices = vertices;
        m.RecalculateNormals();
        m.RecalculateBounds();
    }

    public void CreateParticle(Vector3 position, int idListPos)
    {
        GameObject particle = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        particle.AddComponent<Rigidbody>(); 

        particle.transform.position = position;
        particle.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
        particles.Insert(idListPos, particle);
    }

    public void CreateTriangle(int id1, int id2, int id3)
    {
        Matrix4x4 world2Object = go.transform.worldToLocalMatrix;
        FracTriangle newTriangle = new FracTriangle(world2Object, particles[id1], particles[id2], particles[id3], mu, lambda, phi, psi);
        triangles.Add(newTriangle); 
    }

    public void AttachSpring(int iRb1, int iRb2)
    {
        HingeJoint joint = particles[iRb1].AddComponent<HingeJoint>();
        joint.connectedBody = particles[iRb2].GetComponent<Rigidbody>();

        JointSpring js = joint.spring;
        js.spring = 100;
        js.damper = 10;
        js.targetPosition = 0;

        joint.spring = js;
        joint.useSpring = true;
    }

    public void Fracture(Vector3 force)
    {
        // 1. Find what links get broken

        // 2. duplicate particles where the link breaks (either one of the particles)
        // 3. Re-evaluate links with the new particle
        // 4. apply force to the new particle
        // 5. Recalculate the mesh 
    }

    public void DuplicateParticle(int idP2Copy, int idNewP, ref GameObject newP)
    {
        GameObject iniP = particles[idP2Copy];
        Vector3 pos = iniP.transform.position;
        Vector3 vecPos = new Vector3(pos.x, pos.y, pos.z);

        CreateParticle(vecPos, idNewP);
    }
}
