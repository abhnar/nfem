#include "efem.h"

EFEM::EFEM(){
		epsr[0]=1; epsr[1]=4.3; epsr[2]=2.3;
		mur[0]=1;  mur[1]=1; mur[2]=1;
		tmr =Timer();
	}


	void EFEM::load(string fname){
		
		tmr.tick();

		me.readMesh(fname);
		me.dsplayMeshInfo();

		tmr.tock();

	}

	void EFEM::start(double f) {
		std::cout<<"Portmode Analysis\n";
		tmr.tick();
		ipportmode=portModexy(2, f);
		tmr.tock();

		cout << "Kz\tCut off (GHz)" <<endl<<"-----------------------"<< endl << ipportmode.kz << "\t" << ipportmode.fc<<endl<<endl;
		int cnt = 0;
		for (int i = 0; i < me.tritype.m; i++) {
			if (me.tritype[i] == 2) cnt++;

		}
		cout << cnt << endl;
	}

	/* Extracts single mode in xy plane */

	EFEM::Portmodes EFEM::portModexy(int fno,double f){
		Portmodes pm;
		int edge_seq[3][2]={{0, 1}, {1, 2},{2 ,0}};
		
		Mesh2D port=me.getLocal2DMesh(fno);
	
		int NE=port.nedges;
		int NN= port.nodes.m;

		dmatrix <double> Stt(NE,NE);
		dmatrix <double> Ttt(NE,NE);
		dmatrix <double> Tzt(NN,NE);
		dmatrix <double> Tzz(NN,NN);
		

		double lm,ln;
		double I2=0,I3=0,I4=0,I5=0;
		dmatrix<int> poly;
		poly.zeros(3,4);
		
        f=f*1e9;
        double k0=2*pi/(light_velocity/f);
        cout<<k0<<endl;
        


		double integx=0,integy=0,integx2=0,integy2=0,integxy=0,A,val=0;
		for(int e=0;e<port.triangles.m;e++){
			poly.zeros();
			poly.set(0,1,1);
			double x[3];
			double y[3];
			x[0]=port.nodes(port.triangles(e,0),0);
			x[1]=port.nodes(port.triangles(e,1),0);
			x[2]=port.nodes(port.triangles(e,2),0);
			y[0]=port.nodes(port.triangles(e,0),1);
			y[1]=port.nodes(port.triangles(e,1),1);
			y[2]=port.nodes(port.triangles(e,2),1);
			double a[3]={x[1]*y[2]-y[1]*x[2],x[2]*y[0]-y[2]*x[0],x[0]*y[1]-y[0]*x[1]};
			double b[3]={y[1]-y[2],y[2]-y[0],y[0]-y[1]};
			double c[3]={x[2]-x[1],x[0]-x[2],x[1]-x[0]};
			A=A=abs(x[0]*(y[1]-y[2])+x[1]*(y[2]-y[0])+x[2]*(y[0]-y[1]))/2.0;

			//cout<<port.nodes(port.triangles(e,0),0)<<" "<<port.nodes(port.triangles(e,0),1)<<endl;
			integx= Util::integratePlynomialxy(poly,x[0],x[1],x[2],y[0],y[1],y[2]);
			poly.zeros();
			poly.set(1,1,1);
			integy=Util::integratePlynomialxy(poly,x[0],x[1],x[2],y[0],y[1],y[2]);
			
			poly.zeros();
			poly.set(0,2,1);
			integx2= Util::integratePlynomialxy(poly,x[0],x[1],x[2],y[0],y[1],y[2]);
			
			poly.zeros();
			poly.set(1,2,1);
			integy2= Util::integratePlynomialxy(poly,x[0],x[1],x[2],y[0],y[1],y[2]);

			poly.zeros();
			poly.set(2,0,1);
			integxy= Util::integratePlynomialxy(poly,x[0],x[1],x[2],y[0],y[1],y[2]);
			
			for(int m=0;m<3;m++){
				int i=edge_seq[m][0];
				int j=edge_seq[m][1];
				lm=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
				double Am=a[i]*b[j]-a[j]*b[i];
				double Bm=c[i]*b[j]-c[j]*b[i];
				double Cm=a[i]*c[j]-a[j]*c[i];
				double Dm=-Bm;
				
				for(int n=0;n<3;n++){
					int i=edge_seq[n][0];
					int j=edge_seq[n][1];
					ln=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
					double An=a[i]*b[j]-a[j]*b[i];
					double Bn=c[i]*b[j]-c[j]*b[i];
					double Cn=a[i]*c[j]-a[j]*c[i];
					double Dn=-Bn;
					i=m;j=n;
					I2 = port.signs(e,m)*port.signs(e,n)*((Bm*Bn*lm*ln)/(16*A*A*A*A)-(Bm*Dn*lm*ln)/(16*A*A*A*A) - (Bn*Dm*lm*ln)/(16*A*A*A*A) + 
						(Dm*Dn*lm*ln)/(16*A*A*A*A))*A;

					I3=  port.signs(e,m)*port.signs(e,n)*lm*ln*((Am*An)*A + (Cm*Cn)*A+ (Bm*Bn*integy2)+ (Dm*Dn*integx2) +
					    (Am*Bn*integy) + (An*Bm*integy) + (Cm*Dn*integx) + (Cn*Dm*integx))/(16*A*A*A*A);

					I4=(An*b[i]*port.signs(e,n)*ln)/(8*A*A*A)*A + (Cn*c[i]*port.signs(e,n)*ln)/(8*A*A*A)*A + 
					   (Bn*b[i]*port.signs(e,n)*ln*integy)/(8*A*A*A) + (Dn*c[i]*port.signs(e,n)*ln*integx)/(8*A*A*A);
					 I5=(a[i]*a[j])/(4*A*A)*A + (a[i]*b[j])/(4*A*A)*integx + (a[j]*b[i]*integx)/(4*A*A) + (a[i]*c[j]*integy)/(4*A*A) + (a[j]*c[i]*integy)/(4*A*A) 
					    + (b[i]*b[j]*integx2)/(4*A*A) + (c[i]*c[j]*integy2)/(4*A*A) + (b[i]*c[j]*integxy)/(4*A*A) + (b[j]*c[i]*integxy)/(4*A*A);
					
					 Stt.set(port.tri2edges(e,m),port.tri2edges(e,n),Stt(port.tri2edges(e,m),port.tri2edges(e,n))+I2-k0*k0*epsr[port.tritype[e]]*I3);
					 
					 Ttt.set(port.tri2edges(e,m),port.tri2edges(e,n),Ttt(port.tri2edges(e,m),port.tri2edges(e,n))+I3);

					 Tzt.set(port.triangles(e,i),port.tri2edges(e,n),Tzt(port.triangles(e,i),port.tri2edges(e,n))+I4);

					
					 Tzz.set(port.triangles(e,i),port.triangles(e,j), Tzz(port.triangles(e,i),port.triangles(e,j))+(b[i]*b[j] + c[i]*c[j])/(4*A*A)*A
					 						 - k0*k0*epsr[port.tritype[e]]*I5);
					
					val=val+  I3;


				}
				

			}

		
		}

		
	

	
		
		Stt.deleterows(port.pecedges);
		Stt.deletecols(port.pecedges);
		
		
		Ttt.deleterows(port.pecedges);
		Ttt.deletecols(port.pecedges);

		Tzt.printtofile("Tzt1.txt");

		Tzt.deleterows(port.pecnodes);
		Tzt.deletecols(port.pecedges);
		
		
		Tzz.deleterows(port.pecnodes);
		Tzz.deletecols(port.pecnodes);

		
		//Ttt.printtofile("Ttt.txt");
		//Tzz.printtofile("Tzz.txt");
		//Tzt.printtofile("Tzt.txt");
		//Stt.printtofile("Stt.txt");
		
		

		

		dmatrix <double> Ttz;
		
		Ttz = Tzt.transpose();

		//Ttz.printtofile("Ttz.txt");
		

		/*int info= LAPACKE_dsgesv ( LAPACK_COL_MAJOR, 111, 111, Stt.head, 111,ipiv.head, Ttt.head, 111, x.head, 111, iter);*/
		dmatrix <double> tmp;
		tmp = Tzz.solve(Tzt);


	

		


		/*cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, 
					Ttz.m,tmp.n,Ttz.n,  1, Ttz.head, Ttz.m, tmp.head, tmp.m, 0,c.head,c.m);*/
		dmatrix <double> Btt=Ttz*tmp-Ttt;
	

		//tmp.printtofile("x.txt");
		//Btt.printtofile("Btt.txt");
		dmatrix <double> C= Btt.solve(Stt);
		//C.printtofile("C.txt");
		
		dmatrix<double>  wr(C.n,1);
		dmatrix<double>  ep(C.n, C.n);
	
		dmatrix<double> vr (C.n,C.n);
		dmatrix<double> ev(C.n, 1);

		wr=C.eig(vr);
		cmatrix res;
		res =wr;
		
		cmatrix rr(C.m,1);

		//wr.printtofile("lambda.txt");
		//vr.printtofile("V.txt");
		vzSqrt( res.m, res.head, rr.head);

		//rr.print(0,5);
		wr=rr.imag();
		double smallest_eval=100000;
		int smindex = -1;
		for(int i=0;i<wr.m;i++){
			if(wr(i,0)==0){
				//rr.print(i,i+1);
				if(rr.head[i].real<smallest_eval){
					smallest_eval=rr.head[i].real;
					smindex=i;
				}
			}

		}
		
		//wr.print(90,wr.m);
		if ((k0 * k0 - smallest_eval * smallest_eval) < 0) {
			cout<<"No Propagating Mode"<<endl;
			exit(0);
		}else{
			
			memcpy(ev.head, vr.head + (smindex)*C.n,C.n*sizeof(double));
			//cout << "Norm(ev): " << cblas_dnrm2(C.n, ev.head, 1) << endl;
			//cout <<"Max(ev): "<< ev.max()<< endl;
			//cout << "Min(ev): " << ev.min() << endl;
			//ev.printtofile("V.txt");
			//cout << ev << endl;
			//cout << setprecision(20) << smallest_eval << " " << smindex << endl;
			pm.mode = ev;
			pm.kz = smallest_eval;
			pm.fc = sqrt(k0 * k0 - smallest_eval * smallest_eval)/(2*pi)*light_velocity/1e9;
		}
		
		//cout<<sqrt(k0*k0-smallest_eval*smallest_eval)<<endl;

		return pm;
		
	}
