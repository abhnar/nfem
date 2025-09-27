#pragma once
#include<iostream>
#include<list>
#include<vector>
#include "mkl.h"
#include <iomanip>  

using namespace std;
template <class T>
class dmatrix{
public:

	T *head;

	int m;
	int n;

	dmatrix(){

		head = NULL;

		m=0;
		n=0;

	}
	dmatrix(int mm,int nn){

		head = new T[mm*nn];
		//head=(T *)malloc(mm*nn, sizeof( T ), 64);
		memset(head, 0, sizeof(T) * (nn*mm));

		m=mm;
		n=nn;

	}
	void zeros(int mm,int nn){
		if(head!=NULL)
			free(head);

		head = new T[mm*nn];
		memset(head, 0, sizeof(T) * (nn*mm));

		m=mm;
		n=nn;

	}
	void sequence(){
		if(head!=NULL){
			for(int i=0;i<m*n;i++)
				head[i]=i;
		}else{
			cerr<<"Please initialize the matrix with zeros(int,int)"<<endl;
		}



	}
	void zeros(){
		if(head!=NULL)
			memset(head, 0, sizeof(T) * (n*m));
		else{
			wcerr<<"Please use zeros(int,int)"<<endl;
			exit(1);

		}


	}
	void reset(int mm,int nn){
		if(head!=NULL)
			delete[] (head);
		n=nn;
		m=mm;
		memset(head, 0, sizeof(T) * (n*m));


	}
	

	~dmatrix(){
		delete[] (head);
	}

	const dmatrix &operator=(const dmatrix &a)  {
		m=a.m;
		n=a.n;
		head=new T[m*n];
		for(int i=0;i<m*n;i++){
			head[i]=a.head[i];
		}
		return a;
	}


	/*dmatrix& operator=(dmatrix& a) {
		m = a.m;
		n = a.n;
		head = new T[m * n];
		for (int i = 0; i < m * n; i++) {
			head[i] = a.head[i];
		}
		return a;
	}*/

	inline T &operator()(int index, int y){
		return head[index+y*m]; 
	}
	inline T &operator[](int index){
		return head[index]; 
	}

	inline void set(int x, int y, T val){
		head[x+y*m]=val; 
	}
	inline  T  get(int x, int y){
		return head[x+y*m]; 
	}
	void copyfromvector(const vector<T> &lst){
		head = new T[lst.size()];
		copy(lst.begin(),lst.end(),head);
		m= (int)lst.size();
		n=1;
	}
	void copyfromlist(const list <T> &lst){
		head = new T[lst.size()];
		copy(lst.begin(),lst.end(),head);
		m= static_cast<int>(lst.size());
		n=1;
	}

	void print(){
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				cout<<head[i+m*j]<<" ";
			}
			cout<<endl;
		}
		
	}
	void print(int a,int b){
		for(int i=a;i<b;i++){
			
			cout<<head[i]<<endl;
		}
		
	}
	void printelems(){
		for(int i=0;i<m*n;i++){
			
			cout<<head[i]<<" ";
			
			
		}
		
	}
	
	void printtofile(string fn){
		ofstream myfile;
		myfile.open (fn);
		for(int i=0;i<m;i++){
			for(int j=0;j<n-1;j++){
				myfile<<setprecision(30)<<head[i+j*m]<<",";
			}
			myfile<<setprecision(30)<<head[i+(n-1)*m]<<"\n";
		}
		myfile.close();
		
	}
	friend ostream& operator<<(ostream& stream, const dmatrix& matrix){
		if(matrix.head!=NULL){



			stream<<matrix.m<<"x"<<matrix.n<<" Matrix"<<endl;

			
		}
		else
			stream<<"[]"<<endl;

		return stream;

	}

	void deleterows( dmatrix<int>& rows){
		
		int r1=0, r2, count=0;
		int ix=0;
		int p=0;
		T *ptr= new T[n*(m-rows.m)];
		
		for(int i=0;i<n;i++){
			r1=i*m;
			for(int j=0;j<rows.m;j++){
				
				r2=rows.head[j]+i*m;
				p=r2-r1;
				//cout<<r1<<" "<<r2<<" "<<p<<endl;
				
				memcpy(ptr+count,head+r1,p*sizeof(T));
				r1=r2+1;
				count=count+p;
				
				//cout<<*(head+r1)<<endl;
				
				

			}
			r1=rows.head[rows.m-1]+i*m+1;
			r2=(i+1)*(m);
			p=r2-r1;
			//cout<<r1<<" "<<r2<<" "<<p<<endl;
			
			memcpy(ptr+count,head+r1,p*sizeof(T));
			count=count+p;

			
			
			
			

		}
		//cout<<count<<" "<<n*(m-rows.m)<<endl;
		free(head);
		head=ptr;
		m=m-rows.m;
	}
	void deleterow(int r){
		int el=r;
		int ix=0;
		int p=0;
		for(int i=0;i<n-1;i++){
			el=r+m*i;
			memmove(head+el-i,head+el+1,(n)*sizeof(T));
		}
		int i=n-1;
		el=r+m*i;
		memmove(head+el-i,head+el+1,(m-r)*sizeof(T));

		m=m-1;
		
	}

	void deletecols(const dmatrix<int>& cols){
		int n1=n-cols.m;
		int count=0,ix=0;
		for(int i=0;i<n;i++){
			if(cols.head[count]==i){
				//cout<<cols.head[count]<<endl;
				count++;
				if(count==cols.m-1) count=0;
			}else if(ix!=i){
				memcpy(head+m*ix,head+m*i, m*sizeof(T));
				ix++;
			}else{
				ix++;
			}

		}

		n=n1;
	}
	dmatrix <T> transpose(){
		dmatrix <T> A(n,m);
		
		mkl_domatcopy('C', 'T', m, n, 1, head, m,A.head, A.m);
		return A;
	}

	dmatrix <T>solve(dmatrix<T>& RHS){
		lapack_int * iter;
		iter = new lapack_int;
		dmatrix <lapack_int> ipiv(m,n);
		dmatrix <T> x(RHS.m,RHS.n);
		x.zeros();
		ipiv.zeros();
		int info= LAPACKE_dsgesv ( LAPACK_COL_MAJOR, m, RHS.n, head, m,ipiv.head, RHS.head, RHS.m, x.head, x.m, iter);
		return x;
	}

	dmatrix operator* (const dmatrix& B) const
	{
		dmatrix<double> C(m,B.n);
		cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, 
			m,B.n,n,  1, head, m, B.head, B.m, 0,C.head,C.m);
		return C;
	}
	dmatrix operator- (const dmatrix& B) const
	{
		dmatrix<double> C(m,n);
		mkl_domatadd ('C', 'N', 'N', 
			m ,n, 1, head, m, -1, B.head, B.m, C.head, C.m);
		return C;
	}
	dmatrix operator+ (const dmatrix& B) const
	{
		dmatrix<double> C(m,n);
		mkl_domatadd ('C', 'N', 'N', 
			m ,n, 1, head, m, +1, B.head, B.m, C.head, C.m);
		return C;
	}
	dmatrix eig(dmatrix& B){
		cout << LAPACK_ROW_MAJOR << endl;
		dmatrix<double>  wr(n,1);
		dmatrix<double>  wi(n,1);
		dmatrix<double> vl (n,n);
		dmatrix<double> vr (n,n);
		dmatrix<double> epair(n, n+1);
		LAPACKE_dgeev( LAPACK_COL_MAJOR, 'N', 'V', n, 
			head, m,  wr.head, wi.head, vl.head, vl.m, vr.head, vr.m );
		B=vr;
		return wr;
	}

	T maxx() {
		int maxpos = 0;
		T maximum = head[0];
		for (int i = 1; i < m * n; i++) {
			if (maximum < head[i])
				maximum=head[i];
		}
		return maximum;
	}

	T minn() {
		int maxpos = 0;
		T minimum = head[0];
		for (int i = 1; i < m * n; i++) {
			if (minimum > head[i])
				minimum = head[i];
		}
		return minimum;
	}

	double norm2() {
		int len;
		if (m == 1) {
			len = n;
		}
		else if (n == 1) {
			len = m;
		}
		else {
			cout << "This is a vec norm" << endl;
			return -1;
		}
		double nrm = cblas_dnrm2(len, head, 1);
		return nrm;
	}

	
};

class cmatrix{
public:

	MKL_Complex16 *head;

	int m;
	int n;

	cmatrix(){

		head = NULL;
		m=0;
		n=0;

	}

	cmatrix(int mm,int nn){

		head = new MKL_Complex16[mm*nn];
		//head=(MKL_Complex16 *)malloc(mm*nn, sizeof( MKL_Complex16 ), 64);
		memset(head, 0, sizeof(MKL_Complex16) * (nn*mm));

		m=mm;
		n=nn;

	}
	void zeros(){
		if(head!=NULL)
			memset(head, 0, sizeof(MKL_Complex16) * (n*m));
		else{
			wcerr<<"Please use zeros(int,int)"<<endl;
			exit(1);

		}


	}
	void reset(int mm,int nn){
		if(head!=NULL)
			delete[] (head);
		n=nn;
		m=mm;
		memset(head, 0, sizeof(MKL_Complex16) * (n*m));


	}
	~cmatrix(){
		delete[] (head);
	}

	dmatrix <double> imag() {
		dmatrix <double> a(m,n);
		
		for(int i=0;i<m*n;i++){
			a.head[i]=head[i].imag;
			
		}
		return a;
	}

	dmatrix <double> real() {
		dmatrix <double> a(m,n);
		
		for(int i=0;i<m*n;i++){
			a.head[i]=head[i].real;
			
		}
		return a;
	}

	cmatrix &operator=(cmatrix &a) {
		m=a.m;
		n=a.n;
		head=new MKL_Complex16[m*n];
		for(int i=0;i<m*n;i++){
			head[i].real=a.head[i].real;
			head[i].imag=a.head[i].imag;
		}
		return a;
	}
	dmatrix<double> &operator=(dmatrix<double> &a) {
		m=a.m;
		n=a.n;
		head=new MKL_Complex16[m*n];
		for(int i=0;i<m*n;i++){
			head[i].real=a.head[i];
			head[i].imag=0;
		}
		return a;
	}

	const dmatrix<double> &operator=(const dmatrix<double> &a) {
		m=a.m;
		n=a.n;
		head=new MKL_Complex16[m*n];
		for(int i=0;i<m*n;i++){
			head[i].real=a.head[i];
			head[i].imag=0;
		}
		return a;
	}




	void print(int a,int b){
		for(int i=a;i<b;i++){
			
			cout<<"("<<head[i].real<<","<<head[i].imag<<")"<<endl;
		}
		
	}
	
};