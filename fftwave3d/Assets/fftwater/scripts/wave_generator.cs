using UnityEngine;
using System.Collections;

[ExecuteInEditMode]
[RequireComponent(typeof(MeshFilter))] 
[RequireComponent(typeof(MeshRenderer))]
public class wave_generator : MonoBehaviour {
	public int res_z=8;
	public int res_x=8;
	public float time_scale = 1.0f;
	public float cycle = 200;
	public float gravity = 9.8f;
	public float A =1.0f;
	public Vector2 wind = new Vector2 (1, 0);

	private int real_z = 100;
	private int real_x = 100;
	private float grid_step_x = 0;
	private float grid_step_z = 0;
	private Vector3 [] displace_map;
	public struct Complex{
		public Complex (float real,float image){
			real_=real;
			image_=image;
		}
		public Complex Conj(){
			return new Complex (real_, -image_);
		}
		public static Complex operator *(Complex lhs, Complex rhs){
			return new Complex (lhs.real_ * rhs.real_ - lhs.image_ * rhs.image_,
			                   lhs.image_ * rhs.real_ + lhs.real_ * rhs.image_);
		}
		public static Complex operator +(Complex lhs, Complex rhs){
			return new Complex (lhs.real_ + rhs.real_, lhs.image_ + rhs.image_);
		}
		public override string ToString(){
			return "("+real_.ToString()+","+image_.ToString()+")";
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
	Vector2 calc_k(int m,int n){
		return new Vector2 (Mathf.PI*(2*m-res_x)/real_x,Mathf.PI*(2*n-res_z)/real_z);
	}
	float gauss_random(float mean,float stand_deviation){
		return 1;
	}
	Complex calc_htilde_0(Vector2 k){
		float ei = gauss_random (0, 1);
		float er = gauss_random (0, 1);
		float ph_k = PhillipsTerm (k);
		float temp = Mathf.Sqrt (ph_k / 2);
		return new Complex (temp * er, temp * ei);
	}
	private float PhillipsTerm(Vector2 k){
		float k_mag = k.magnitude;
		float L = Vector2.Dot(wind , wind) / gravity;
		float term1 = (Mathf.Exp (-1.0f / Mathf.Pow (k_mag * L, 2))) / (Mathf.Pow (k_mag, 4));
		float term2 = Mathf.Pow (Vector2.Dot (k.normalized, wind.normalized), 2);
		return A*term1*term2;
	}
	//exp(ix) = cosx + isinx
	public static Complex Expi(float x){
		return new Complex (Mathf.Cos (x), Mathf.Sin (x));
	}  
	//wk=sqrt(gk)
	//[[wk/w0]]*w0
	//w0=2*PI/T
	public float calc_dispersion(Vector2 k){
		float k_mag = k.magnitude;
		float omega_0 = 2 * Mathf.PI / cycle;
		return (Mathf.Floor(Mathf.Sqrt(gravity*k_mag)/omega_0))*omega_0;
	}
	// h0(k)*exp(i*wk*t)
	//+h0(-k).conj()*exp(-i*wk*t)
	Complex calc_htilde(int m,int n,float t){
		Vector2 k = calc_k (m, n);
		float wt = calc_dispersion (k) * t;
		return calc_htilde_0 (k) * Expi (wt) + calc_htilde_0 (-k).Conj () * Expi (-wt);
	}
	void dft(Complex [] data_in,Complex [] data_out,int start,int stride){
		int N = (data_in.Length-start)/stride;
		Complex [] temp_res = new Complex[N];
		for (int k=0; k<N; k++) {
			temp_res[k] = new Complex(0,0);
			for (int j=0;j<N;j++){
				temp_res[k]=temp_res[k]+Expi(2*Mathf.PI*j*k/N)*data_in[start+j*stride];
			}
			data_out[start+k*stride]=temp_res[k];
		}
	}
	//dft
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
			//fft (h_tilde,h_tilde,m*res_z,1);
			dft (h_tilde,h_tilde,m*res_z,1);
		}
		for (int n=0; n<res_z; n++) {
			//fft (h_tilde,h_tilde,n,res_z);
			dft (h_tilde,h_tilde,n,res_z);
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
