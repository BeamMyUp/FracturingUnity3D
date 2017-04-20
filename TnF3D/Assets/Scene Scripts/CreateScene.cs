using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CreateScene : MonoBehaviour {

    List<Particles_Mesh> objects = new List<Particles_Mesh>(); 
    public int meshPoints;
    
    void CreateSquareMesh(ref GameObject go, ref List<GameObject> particles)
    {
        // Create a basic squared mesh
        Vector3[] vertices = new Vector3[particles.Count];
        
        for (int i = 0; i < particles.Count; ++i)
        {
            vertices[i] = particles[i].transform.position; 
        }

        List<int> indices = new List<int>();
        for (int y = 0; y < meshPoints - 1; ++y)
            for (int x = 0; x < meshPoints - 1; ++x)
            {
                int val_ul = x + (y * meshPoints);
                int val_ur = x + (y * meshPoints) + 1;
                int val_bl = x + ((y + 1) * meshPoints);
                int val_br = x + ((y + 1) * meshPoints) + 1;

                indices.Add(val_ul);
                indices.Add(val_bl);
                indices.Add(val_br);

                indices.Add(val_ul);
                indices.Add(val_br);
                indices.Add(val_ur);
            }

        // Create the mesh
        Mesh msh = new Mesh();
        msh.vertices = vertices;
        msh.triangles = indices.ToArray();
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

    void CreateSpringGrid(ref List<GameObject> particles)
    {
        int lowerlim = meshPoints / 2;
        int higherlim = lowerlim + (meshPoints % 2); 

        // create all spheres
        for (int y = -lowerlim; y < higherlim; ++y)
            for (int x = -lowerlim; x < higherlim; ++x)
            {
                GameObject sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                sphere.AddComponent<Rigidbody>();

                sphere.transform.position = new Vector3(x, y, 0);
                sphere.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
                sphere.AddComponent<DragAndDrop>(); 
                particles.Add(sphere); 
            }

        // Assign springs
        for (int j = 0; j < particles.Count; ++j)
        {
            bool inLastCol = (j + 1) % meshPoints == 0;
            bool inLastRow = j >= meshPoints * (meshPoints - 1);

            // horizonal springs
            if (!inLastCol) attachSpring(ref particles, j, j + 1); 

            // vertical springs
            if (!inLastRow) attachSpring(ref particles, j, j + meshPoints); 

            // diagonal springs
            if(!inLastCol && !inLastRow) attachSpring(ref particles, j, j + meshPoints + 1);
        }
    }

    void attachSpring(ref List<GameObject> particles, int iRb1, int iRb2)
    {
        HingeJoint joint = particles[iRb1].AddComponent<HingeJoint>();
        joint.connectedBody = particles[iRb2].GetComponent<Rigidbody>();
        
        JointSpring js = joint.spring;
        js.spring = 10;
        js.damper = 3;
        js.targetPosition = 0;
        
        joint.spring = js;
        joint.useSpring = true;
    }

	// Use this for initialization
    void Awake()
    {
        // Create floor 
        GameObject floorPlane = GameObject.CreatePrimitive(PrimitiveType.Plane);
        floorPlane.transform.position = new Vector3(0, -meshPoints/2 - 1, 0);
        floorPlane.transform.localScale = new Vector3(2, 1, 2);

        // Create a first Particles_Mesh
        List<GameObject> particles = new List<GameObject>();
        GameObject mesh = new GameObject("Tearable Object"); 

        CreateSpringGrid(ref particles);
        CreateSquareMesh(ref mesh, ref particles);

        Particles_Mesh initObject = new Particles_Mesh(ref particles, ref mesh);
        objects.Add(initObject);
    }
    
	// Update is called once per frame
	void Update ()
    {
        for(int i = 0; i < objects.Count; ++i)
        {
            objects[i].Update(); 
        }
	}
}
