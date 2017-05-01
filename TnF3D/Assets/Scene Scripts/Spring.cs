using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Spring : MonoBehaviour {

    Rigidbody p1 = null;
    Rigidbody p2 = null;

    // Spring Stiffness
    public static double k = 700;

    // Spring damping
    public static double c = 10;

    // Rest Length
    double l0 = 0; 

    public Rigidbody P1
    {
        get { return p1; }
    }

    public Rigidbody P2
    {
        get { return p2; }
    }

    public void Initialize(Rigidbody P1, Rigidbody P2)
    {
        if(p1 != null)
        {
            p1.gameObject.GetComponent<FracParticle>().Springs.Remove(this);
        }
        if(p2 != null)
        {
            p2.gameObject.GetComponent<FracParticle>().Springs.Remove(this);
        }

        p1 = P1;
        p2 = P2;
        l0 = Vector3.Distance(p1.transform.position, p2.transform.position);

        p1.gameObject.GetComponent<FracParticle>().Springs.Add(this);
        p2.gameObject.GetComponent<FracParticle>().Springs.Add(this);
    }
    
    public GameObject GetOtherPoint(GameObject p)
    {
        if (p == p1.gameObject)
            return p2.gameObject;
        if (p == p2.gameObject)
            return p1.gameObject;

        return null; 
    }

    void CalculateSpringForce()
    {
        Vector3 l = p1.transform.position - p2.transform.position;
        Vector3 l_ = p1.velocity - p2.velocity;
        double l1 = l.magnitude;

        double fa = -((k * (l1 - l0)) + (c * (Vector3.Dot(l, l_) / l1))) / l1;

        Vector3 forceA = l * (float)fa;
        Vector3 forceB = -forceA;

        p1.AddForce(forceA);
        p2.AddForce(forceB);
    }

    // HACK: This is a test function, not to use really.
    public bool IsBreaking(Vector3 planeN, GameObject point)
    { 
        if (l0 > l0 + 0.2f)
        {
            planeN = p2.position - p1.position;
            if (p1.velocity.magnitude > p2.velocity.magnitude)
            {
                point = p2.gameObject;
            }
            else
            {
                point = p1.gameObject;
            }

            return true; 
        }

        return false; 
    }

    // Update is called once per frame
    void Update () {
        CalculateSpringForce();
	}
}
