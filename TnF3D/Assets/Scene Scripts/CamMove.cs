using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CamMove : MonoBehaviour {

    float transSpeed = 0.05f;
    float rotSpeed = 0.01f;
    Vector2 initMousePos; 

	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
        // rotate
		if(Input.GetKey(KeyCode.LeftAlt)){
            if(Input.GetMouseButtonDown(0)){
                initMousePos = Input.mousePosition; 
            }

            if(Input.GetMouseButton(0)){
                Vector2 newMousePos = Input.mousePosition;
                Vector2 theta = newMousePos - initMousePos;

                transform.Rotate(Vector3.up*theta.x*rotSpeed, Space.World);
                transform.Rotate(Vector3.right*-theta.y*rotSpeed); 
            }
        }
        Vector3 translation = new Vector3(); 

        // translate left
        if (Input.GetKey(KeyCode.LeftArrow)){
            translation -= transform.right * transSpeed;  
        }
        // translate right
        if (Input.GetKey(KeyCode.RightArrow)){
            translation += transform.right * transSpeed; 
        }
        // move forward
        if (Input.GetKey(KeyCode.UpArrow)){
            translation += transform.forward * transSpeed;
        }
        // move backward
        if (Input.GetKey(KeyCode.DownArrow)){
            translation -= transform.forward * transSpeed; 
        }

        transform.Translate(translation); 

	}
}
