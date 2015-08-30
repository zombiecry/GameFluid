using UnityEngine;
using System.Collections;

[ExecuteInEditMode]
[RequireComponent(typeof(MeshFilter))] 
[RequireComponent(typeof(MeshRenderer))]
public class wave_generator : MonoBehaviour {
	public int res_z=8;
	public int res_x=8;
	public float time_scale = 1.0f;
	private int real_z = 100;
	private int real_x = 100;
	private float grid_step_x = 0;
	private float grid_step_z = 0;
	private Vector3 [] displace_map;
	struct Complex{
		public Complex (float real,float image){
			real_=real;
			image_=image;
		}
		public float real_;
		public float image_;
	};
	//fft related
	private Complex [] h_tilde;
	// Use this for initialization
	void Start () {
		//init
		displace_map = new Vector3[res_x * res_z];
		for (int i=0; i<displace_map.Length; i++) {
			displace_map[i]= new Vector3 (0,0,0);
		}
		h_tilde = new Complex[res_x * res_z];
		//-----------------------------------------
		//init the mesh
		Vector3 [] verts = new Vector3[(res_x + 1) * (res_z + 1)];
		int     [] tris  = new int[res_x*res_z*6];
		grid_step_x = (float)real_x / (float)res_x;
		grid_step_z = (float)real_z / (float)res_z;
		for (int i=0;i<=res_x;i++) for (int j=0;j<=res_z;j++){
			int vert_idx = i*(res_z+1)+j;
			verts[vert_idx] = new Vector3(i*grid_step_x,0,j*grid_step_z);
			if (i>=res_x || j>=res_z) continue;
			int grid_idx = i*(res_z)+j;
			tris[grid_idx*6 + 0] = vert_idx;
			tris[grid_idx*6 + 1] = vert_idx+(res_z+1)+1;
			tris[grid_idx*6 + 2] = vert_idx+(res_z+1);
			tris[grid_idx*6 + 3] = vert_idx;
			tris[grid_idx*6 + 4] = vert_idx+1;
			tris[grid_idx*6 + 5] = vert_idx+(res_z+1)+1;
		}
		MeshFilter mf = GetComponent<MeshFilter> ();
		Mesh mesh = new Mesh();
		mesh.vertices = verts;
		mesh.triangles = tris;
		mf.mesh = mesh;
		//-------------------------------------------
	}
	Complex calc_htilde(int m,int n,float t){
		// to do...
		return new Complex (1, 0);
	}

	void fft(Complex [] data_in,Complex [] data_out,int start,int stride){
		for (int i=start; i<data_in.Length; i+=stride) {
			data_out[i]=data_in[i];
		}
	}

	void gen_dis_map(float t){
		for (int m=0;m<res_x;m++) for (int n=0;n<res_z;n++){
			h_tilde[m*res_z+n] = calc_htilde(m,n,t);
		}
		for (int m=0;m<res_x;m++){
			fft (h_tilde,h_tilde,m*res_z,1);
		}
		for (int n=0; n<res_z; n++) {
			fft (h_tilde,h_tilde,n,res_z);
		}
		float sign = 1;
		for (int i=0; i<res_x; i++) for (int j=0; j<res_z; j++) {
			if ((i+j)%2==1) sign = -1;
			displace_map[i*res_z+j] = new Vector3(0,sign * h_tilde[i*res_z+j].real_);
		}
	}
	void update_mesh(){
#if UNITY_EDITOR
		Mesh mesh = GetComponent<MeshFilter>().sharedMesh;
#else
		Mesh mesh =GetComponent<MeshFilter> ().mesh;
#endif
		Vector3 [] verts = mesh.vertices;
		for (int i=0; i<=res_x; i++) for(int j=0;j<res_z;j++){
			int vert_idx = i*(res_z+1)+j;
			verts[vert_idx]=new Vector3(i*grid_step_x,0,j*grid_step_z); 
			verts[vert_idx]+=displace_map[(i%res_x)*res_z+j%res_z];
		}
	}
	// Update is called once per frame
	void Update () {
		gen_dis_map (Time.time*time_scale);
		update_mesh ();
	}
}
