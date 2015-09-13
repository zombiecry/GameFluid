#define NO_FT
//#define USE_FFT
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
	public Vector2 wind = new Vector2 (3, 0);

	private int real_z = 8;
	private int real_x = 8;


	private float grid_step_x = 0;
	private float grid_step_z = 0;
	private Vector3 [] displace_map;

	private Complex [] htilde_0;
	private Complex [] htilde_0_1;
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
		public static Complex operator *(Complex lhs, float x){
			return new Complex (lhs.real_ * x,lhs.image_* x);
		}
		public static Complex operator +(Complex lhs, Complex rhs){
			return new Complex (lhs.real_ + rhs.real_, lhs.image_ + rhs.image_);
		}
		public static Complex operator -(Complex lhs, Complex rhs){
			return new Complex (lhs.real_ - rhs.real_, lhs.image_ - rhs.image_);
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
		htilde_0 = new Complex[res_x * res_z];
		htilde_0_1 = new Complex[res_x * res_z];
		for (int m=0;m<res_x;m++) for (int n=0;n<res_z;n++){
			Vector2 k = calc_k (m, n);
			htilde_0[m*res_z+n] = calc_htilde_0(k);
			htilde_0_1[m*res_z+n] = calc_htilde_0 (-k).Conj ();
			//htilde_0[m*res_z+n] = new Complex(1,0);
			
			//htilde_0_1[m*res_z+n] = new Complex(1,0);
		}

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
	Complex gauss_random(){
		float x1, x2, w; 
		do {
			x1 = Random.value;
			x2 = Random.value;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0f );
		w = Mathf.Sqrt((-2.0f * Mathf.Log(w)) / w);
		return new Complex(x1 * w, x2 * w);
	}
	Complex calc_htilde_0(Vector2 k){
		Complex r = gauss_random ();
		float ph_k = PhillipsTerm (k);
		float temp = Mathf.Sqrt (ph_k / 2.0f);
		return r*temp;
	}
	private float PhillipsTerm(Vector2 k){
		float k_mag = k.magnitude;
		if (k_mag < 0.00001) return 0;
		float k_mag2 = k_mag * k_mag;
		float k_mag4 = k_mag2 * k_mag2;

		float k_dot_w = Vector2.Dot (k.normalized, wind.normalized);
		float k_dot_w2 = k_dot_w * k_dot_w;

		float l = wind.magnitude * wind.magnitude;
		float l2 = l * l;

		float damping = 0.001f;
		float l2d = l2 * damping*damping;
		return A * Mathf.Exp (-1.0f / (k_mag2 * l2)) / k_mag4 * k_dot_w2 * Mathf.Exp (-k_mag2 * l2d);
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
		float omega_0 = 2.0f * Mathf.PI / cycle;
		return (Mathf.Floor(Mathf.Sqrt(gravity*k_mag)/omega_0))*omega_0;
	}
	// h0(k)*exp(i*wk*t)
	//+h0(-k).conj()*exp(-i*wk*t)
	void dft(Complex [] data_in,Complex [] data_out,int start,int stride,int N){
		Complex [] temp_res = new Complex[N];
		for (int k=0; k<N; k++) {
			temp_res[k] = new Complex(0,0);
			for (int j=0;j<N;j++){
				temp_res[k]=temp_res[k]+Expi(2*Mathf.PI*j*k/N)*data_in[start+j*stride];
			}
		}
		for (int k=0; k<N; k++) {
			data_out[start+k*stride]=temp_res[k];
		}
	}
	private int fast_log_2(int x){
		int count = 0;
		while ((x/=2)>0)	count++;
		return count;
	}
	private int flip(int x,int n){
		int l = 1;
		int h = n;
		int y = x;
		while (l<h) {
			int bh = x&h;
			int bl = x&l;
			if ((bh>0) != (bl>0)){
				y&=~h;
				y&=~l;
				y|=(~bh)&h;
				y|=(~bl)&l;
			}
			h=h>>1;
			l=l<<1;
		}
		return y;
	}
	static void Swap(ref Complex lhs, ref Complex rhs)
	{
		Complex temp;
		temp = lhs;
		lhs = rhs;
		rhs = temp;
	}

	private void reorder(Complex [] data,int start,int stride,int N){
		//string temp = "";
		for (int i=0; i<N; i++) { 
			int j = flip(i,N>>1);
			if (i<j) {
				/*
				Debug.Log("--------------------------------");
				Debug.Log(data[start+i*stride].ToString());
				Debug.Log(data[start+j*stride].ToString());

				Debug.Log(data[start+i*stride].ToString());
				Debug.Log(data[start+j*stride].ToString());
				Debug.Log("--------------------------------");
				*/
				/*
				temp +="<";
				temp +=i.ToString();
				temp +=":";
				temp+=j.ToString();
				temp+="> ";
				*/
				Swap(ref data[start+i*stride],ref data[start+j*stride]);
			}
		}
		//Debug.Log (temp);
	}
	//dft
	void fft(Complex [] data_in,Complex [] data_out,int start,int stride,int N){
		for (int i=0; i<N; i++) {
			data_out[start+i*stride] = data_in[start+i*stride];
		}
		reorder (data_out,start,stride,N);
		int log_n =fast_log_2(N);
		for (int m=0;m<log_n;m++){					//stage
			int nm=1<<(m+1);
			Complex w = new Complex(1,0);
			Complex wm = Expi(-2*Mathf.PI/(float)nm);
			int freq_num = nm>>1;
			for (int k=0;k<freq_num;k++){			//every freq in cur stage
				for (int i=k;i<N;i+=nm){			//every butterfly in cur freq,start with offset of k
					int p=i;
					int q=i+freq_num;				//freq_num is Nm/2
					int p_pos=start+p*stride;
					int q_pos=start+q*stride;
					Complex t = w*data_out[q_pos];
					data_out[q_pos]=data_out[p_pos]-t;
					data_out[p_pos]=data_out[p_pos]+t;
				}
				w=w*wm;
			}
		}
	}
	Complex calc_htilde(int m,int n,float t){
		Vector2 k = calc_k (m, n);
		float wt = calc_dispersion (k) * t;
		return htilde_0[m*res_z+n] * Expi (wt) + htilde_0_1[m*res_z+n] * Expi (-wt);
	}

	void gen_dis_map(float t){
		for (int m=0;m<res_x;m++) for (int n=0;n<res_z;n++){
			h_tilde[m*res_z+n] = calc_htilde(m,n,t);
			//temp +="<"+m.ToString()+","+n.ToString ()+">";
			//temp +=h_tilde[m*res_z+n].ToString();
			//temp+= calc_k (m, n).ToString();
		}

#if NO_FT
		Complex [] temp1 = new Complex[res_x * res_z];
		for (int x=0; x<res_x; x++) for (int z=0; z<res_z; z++) {
				temp1 [x * res_z + z] = new Complex (0, 0);
				Vector2 vx = new Vector2((x) * real_x/(float)res_x, (z)*real_z/(float)res_z);
				for (int m=0; m<res_x; m++) for (int n=0; n<res_z; n++) {
						Vector2 k = calc_k (m, n);
						Complex b = Expi (Vector2.Dot (k, vx));
						Complex a = h_tilde[m*res_z+n];
						temp1[x*res_z +z]=temp1[x*res_z+z]+a*b;
				}
				displace_map[x*res_z+z] = new Vector3(0,temp1[x*res_z+z].real_*0.03f,0);
		}
#else

		for (int m=0;m<res_x;m++){
#if USE_FFT
			fft (h_tilde,h_tilde,m*res_z,1,res_z);
#else
			dft (h_tilde,h_tilde,m*res_z,1,res_z);
#endif
		}
		for (int q=0; q<res_z; q++) {
#if USE_FFT
			fft (h_tilde,h_tilde,q,res_z,res_x);
#else
			dft (h_tilde,h_tilde,q,res_z,res_x);
#endif
		}
		float sign = 1;
		for (int i=0; i<res_x; i++) for (int j=0; j<res_z; j++) {
			if ((i+j)%2==0) sign=-1;
			else sign=1;
			displace_map[i*res_z+j] = new Vector3(0,sign * h_tilde[i*res_z+j].real_*0.03f,0);
		}
#endif

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
		mesh.vertices = verts;
	}
	// Update is called once per frame
	void Update () {
		gen_dis_map (Time.time*time_scale);
		update_mesh ();
	}
}
