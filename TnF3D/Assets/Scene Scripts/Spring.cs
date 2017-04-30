using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Spring : MonoBehaviour {

    Rigidbody p1 = null;
    Rigidbody p2 = null;

    // Spring Stiffness
    public static double k = 800;

    // Spring damping
    public static double c = 10;

    // Rest Length
    double l0 = 0; 

    public void Initialize(Rigidbody p1, Rigidbody p2)
    {
        this.p1 = p1;
        this.p2 = p2;
        l0 = Vector3.Distance(p1.transform.position, p2.transform.position);
    }


	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
        Vector3 l = p1.transform.position - p2.transform.position;
        Vector3 l_ = p1.velocity - p2.velocity;
        double l1 = l.magnitude;
        
        double fa = -((k * (l1 - l0)) + (c * (Vector3.Dot(l, l_) / l1))) / l1;

        Vector3 forceA = l * (float)fa;
        Vector3 forceB = -forceA;

        p1.AddForce(forceA);
        p2.AddForce(forceB); 
	}
}
