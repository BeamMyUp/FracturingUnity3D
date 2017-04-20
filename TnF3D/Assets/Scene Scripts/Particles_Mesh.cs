using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particles_Mesh
{
    Vector3[] vertices;

    GameObject go;
    List<GameObject> particles = new List<GameObject>();

    public Particles_Mesh(ref List<GameObject> aParticles, ref GameObject aGO)
    {
        particles = aParticles;
        go = aGO;
        vertices = go.GetComponent<MeshFilter>().mesh.vertices; 
    }
	
	public void Update()
    {
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position;
        }

        Mesh m = go.GetComponent<MeshFilter>().mesh;
        m.vertices = vertices;
        m.RecalculateNormals();
        m.RecalculateBounds();
    }

    public void IntraSeparate(Vector3 force)
    {
        // 1. Find what links get broken

        // 2. duplicate particles where the link breaks (either one of the particles)
        // 3. Re-evaluate links with the new particle
        // 4. apply force to the new particle
        // 5. Recalculate the mesh 
    }

    public void DuplicateParticle(int idP, ref GameObject newP)
    {
        GameObject iniP = particles[idP]; 

        newP = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        newP.AddComponent<Rigidbody>();

        Vector3 pos = iniP.transform.position;

        newP.transform.position = new Vector3(pos.x, pos.y, pos.z);
        newP.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
        newP.AddComponent<DragAndDrop>();
    }
}
