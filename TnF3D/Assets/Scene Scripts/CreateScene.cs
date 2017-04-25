using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CreateScene : MonoBehaviour {

    List<GameObject> fracMeshes = new List<GameObject>(); 
    public int meshPoints;

    // TODO: should be refactored. Indices creation shouldn't be independent of vertices/particles creation
    void CreateSquareMesh(FracMesh fm)
    {
        List<int> indices = new List<int>();

        for (int z = 0; z < meshPoints - 1; ++z)
            for (int x = 0; x < meshPoints - 1; ++x)
            {
                int val_ul = x + (z * meshPoints);
                int val_ur = x + (z * meshPoints) + 1;
                int val_bl = x + ((z + 1) * meshPoints);
                int val_br = x + ((z + 1) * meshPoints) + 1;

                indices.Add(val_ul);
                indices.Add(val_bl);
                indices.Add(val_br);
                fm.CreateTriangle(val_ul, val_bl, val_br);

                indices.Add(val_ul);
                indices.Add(val_br);
                indices.Add(val_ur);
                fm.CreateTriangle(val_ul, val_br, val_ur);
            }

        fm.CreateMesh(indices.ToArray()); 
    }

    void CreateSpringGrid(FracMesh fm)
    {
        int lowerlim = meshPoints / 2;
        int higherlim = lowerlim + (meshPoints % 2); 

        // create all spheres
        for (int z = -lowerlim; z < higherlim; ++z)
            for (int x = -lowerlim; x < higherlim; ++x)
            {
                Vector3 pos = new Vector3(x, 0, z);
                fm.CreateParticle(pos, fm.nParticles); 
            }

        // Assign springs
        for (int j = 0; j < fm.nParticles; ++j)
        {
            bool inLastCol = (j + 1) % meshPoints == 0;
            bool inLastRow = j >= meshPoints * (meshPoints - 1);

            // horizonal springs
            if (!inLastCol) fm.AttachSpring(j, j + 1); 

            // vertical springs
            if (!inLastRow) fm.AttachSpring(j, j + meshPoints); 

            // diagonal springs
            if(!inLastCol && !inLastRow) fm.AttachSpring(j, j + meshPoints + 1);
        }
    }

	// Use this for initialization
    void Awake()
    {
        // Create floor 
        GameObject floorPlane = GameObject.CreatePrimitive(PrimitiveType.Plane);
        floorPlane.transform.position = new Vector3(0, -meshPoints/2 - 1, 0);
        floorPlane.transform.localScale = new Vector3(2, 1, 2);

        // Create a first Fracturable Mesh
        GameObject initObject = new GameObject();
        FracMesh fm = initObject.AddComponent<FracMesh>();
        fm.InitializeFracMesh(2650000, 3970000, 264, 397);

        CreateSpringGrid(fm);
        CreateSquareMesh(fm);

        fracMeshes.Add(initObject);
    }
    
	// Update is called once per frame
	void Update ()
    {
        
	}
}
