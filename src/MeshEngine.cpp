#include "meshengine.h"

MeshEngine::MeshEngine(){
        meshfilename="";
        domcount=0;
        boundcount=10;
        
       
       // nodes=NULL;
        tettype=NULL;
        
        signs_tet=NULL;
        signs_tri=NULL;

        
    }


    MeshEngine::~MeshEngine(){
       
        
       
        if(tettype!=NULL){
            free(tettype);
        }
       
        if(tetoftris!=NULL){
            free(tetoftris);
        }
        if(signs_tet!=NULL){
            for(int i=0;i<ntets;i++)
                free(signs_tet[i]);
            free(signs_tet);
        }
        if(signs_tri!=NULL){
            for(int i=0;i<ntris;i++)
                free(signs_tri[i]);
            free(signs_tri);
        }

    }



     Mesh2D  MeshEngine::getLocal2DMesh(int fno){

        dmatrix <int>tri =getTriangles(fno);

        dmatrix <int>pec =getTriangles(1);

        list <int> lst_n;
        list <int> lst_e;
        list <int> lst_pn;
        list <int> lst_pe;
        list <int> lst_pec;
        std::cout<<tri.m<<std::endl;
        for(int i=0;i<tri.m;i++){
            std::cout<<tri[i]<<std::endl;
            for(int j=0;j<3;j++){
                lst_n.push_back(triangles(tri[i],j));
                lst_e.push_back(tri2edge(tri[i],j));
            }
            //cout<<endl;
        }
        

        
        lst_n.sort();
        lst_n.unique();
        lst_e.sort();
        lst_e.unique();

        dmatrix<int> nds;
        nds.copyfromlist(lst_n);
        dmatrix<int> edgs;
        edgs.copyfromlist(lst_e);


        for(int i=0;i<pec.m;i++){
            for(int j=0;j<3;j++){
                lst_pn.push_back(triangles(pec[i],j));
                lst_pe.push_back(tri2edge(pec[i],j));
            }

        }


        lst_pe.sort();
        lst_pe.unique();
        lst_pn.sort();
        lst_pn.unique();

        
        set_intersection(lst_e.begin(), lst_e.end(), lst_pe.begin(), lst_pe.end(),
                          std::back_inserter(lst_pec));
        dmatrix<int> pec_edgs;
        pec_edgs.copyfromlist(lst_pec);

        lst_pec.clear();
        set_intersection(lst_n.begin(), lst_n.end(), lst_pn.begin(), lst_pn.end(),
                          std::back_inserter(lst_pec));
        

        dmatrix<int> pec_nds;
        pec_nds.copyfromlist(lst_pec);

        map<int, int> nodedict;
        map<int, int> edgedict;
        //cout<<nds.m<<endl;
        dmatrix<double> nodes_local(nds.m,3);
        for(int i=0;i<nds.m;i++){
            nodedict.insert(pair <int,int> (nds[i],i));
            
            nodes_local.set(i,0,nodes(nds[i],0));
            nodes_local.set(i,1,nodes(nds[i],1));
            nodes_local.set(i,2,nodes(nds[i],2));
        }

        //cout<<nodedict.size()<<endl;
        dmatrix<int> trinagles_local(tri.m,3);
        dmatrix <int>tritype_local(tri.m,1);
        dmatrix <int>pec_nds_local(pec_nds.m,1);
        dmatrix <int>pec_edgs_local(pec_edgs.m,1);
        
        for(int i=0;i<tri.m;i++){

            trinagles_local.set(i,0,nodedict[triangles(tri[i],0)]);
            trinagles_local.set(i,1,nodedict[triangles(tri[i],1)]);
            trinagles_local.set(i,2,nodedict[triangles(tri[i],2)]);

            tritype_local.set(i,0,tettype[tetoftris[tri[i]]]);

        }
        for(int i=0;i<pec_nds.m;i++){
            pec_nds_local.set(i,0,nodedict[pec_nds[i]]);
        }

        for(int i=0;i<edgs.m;i++){
            edgedict.insert(pair <int,int> (edgs[i],i));
        }
        dmatrix<int> tri2edges_local(tri.m,3);
        for(int i=0;i<tri.m;i++){

            tri2edges_local.set(i,0,edgedict[tri2edge(tri[i],0)]);
            tri2edges_local.set(i,1,edgedict[tri2edge(tri[i],1)]);
            tri2edges_local.set(i,2,edgedict[tri2edge(tri[i],2)]);

        }
        for(int i=0;i<pec_edgs.m;i++){
            pec_edgs_local.set(i,0,edgedict[pec_edgs[i]]);
        }

        int FEdge_seq[3][2]={{0, 1}, {1, 2},{2 ,0}};
        dmatrix<int> signs_local(tri2edges_local.m,3);
        
        for(int i=0;i<tri2edges_local.m;i++){
            for(int j=0;j<3;j++){
                signs_local.set(i,j,(trinagles_local(i,FEdge_seq[j][0])-trinagles_local(i,FEdge_seq[j][1]))<0?-1:1);
            }
        }

        
        Mesh2D facemesh;
        facemesh.nodes=nodes_local;
        facemesh.pecnodes=pec_nds_local;
        facemesh.pecedges=pec_edgs_local;
        facemesh.triangles = trinagles_local;
        facemesh.tri2edges=tri2edges_local;
        facemesh.signs=signs_local;
        facemesh.nedges=edgs.m;
        facemesh.tritype=tritype_local;
        
        
        
        
        return facemesh;
        


    }

     void MeshEngine::displayTagInfo(){


        cout<<endl<<"********************************************"<<endl<<"BOUNDARY INFORMATION"<<endl<<"********************************************"<<endl;
        for(int i=0;i<Boundaries.size();i++){
            cout<<"Boundary Name: "<<Boundaries[i].name<<endl<<"tagId : "<<Boundaries[i].tagId<<endl<<"tag(s): ";
            for(int j=0;j<Boundaries[i].tags.size();j++)
                cout<<Boundaries[i].tags[j]<<" ";
            cout<<endl<<"--------------------------------------------"<<endl;
        }

        cout<<endl<<"############################################"<<endl<<"DOMAIN INFORMATION"<<endl<<"############################################"<<endl;
        for(int i=0;i<Domains.size();i++){
            cout<<"Domain Name: "<<Domains[i].name<<endl<<"tagId : "<<Domains[i].tagId<<endl<<"tag(s): ";
            for(int j=0;j<Domains[i].tags.size();j++)
                cout<<Domains[i].tags[j]<<" ";
            cout<<endl<<"--------------------------------------------"<<endl;
        }
        

    }


    void MeshEngine::dsplayMeshInfo(){
        
        cout<<"#Nodes\t#Edges\t#Tris \t#Tets"<<"\n-------------------------------"<<endl;
        cout<<nnodes<<"\t"<<nedges<<"\t"<<ntris<<"\t"<<ntets<<"\n-------------------------------"<<endl;

        cout<<"Boundaries :";
        for(int i=0;i<Boundaries.size();i++){
            cout<< Boundaries[i].tagId <<" - "<<Boundaries[i].name<<"  ||  ";
        }
        cout<<endl<<"Domains    :";
        for(int i=0;i<Domains.size();i++){
            cout<< Domains[i].tagId <<" - "<<Domains[i].name<<"  ||  ";
        }
        cout<<endl;
        
        
    }
    void  MeshEngine::readMesh(string str){
        meshfilename=str;
        cout << "Getting tag Info"<<endl;
        getTagInfo();
        cout<<"Mesh reading stage - 1";
        getMeshInfo();
        cout<<"\b2";
        getEdgeData();
        cout<<"\b3";
        getTetofTri();
        cout<<"\b4 (Finished)"<<endl;
        
    }


    dmatrix<int> MeshEngine::getTriangles(int t){
        dmatrix <int> tri;
        vector <int> v;
        for(int i=0;i<ntris;i++){
            if(tritype[i]==t)
                v.push_back(i);
        }
        tri.copyfromvector(v);
        return(tri);

    }
    dmatrix<int> MeshEngine::getTetrahedrons(int t){
        dmatrix <int> tet;
        vector <int> v;
        for(int i=0;i<ntets;i++){
            if(tettype[i]==t)
                v.push_back(i);
        }
        tet.copyfromvector(v);
        return(tet);

    }


    void MeshEngine::getVolume() {
        for (int i = 0; i < ntets;  i++) {

        }
        cout << ntets << endl;

    }


    int MeshEngine::getTetofTri(){
        tetoftris= new int[ntris];
        int count=0;
        for(int ix=0;ix<ntris;ix++){
            tetoftris[ix]=-1;
            for(int i=0;i<ntets;i++){
                count=0;
                for(int j=0;j<3;j++){
                    for(int k=0;k<4;k++){
                        if(tetrahedrons(i,k) == triangles(ix,j)){
                            count++;
                            break;
                        }
                    }

                }
                if(count==3){
                    tetoftris[ix]=i;
                    break;
                }

            }

        }
        return 0;

    }
    void MeshEngine::signs(){
        int Edge_seq[6][2]={{0, 1}, {0, 2},{0 ,3}, {1, 2}, {3 ,1}, {2, 3}};
       
        signs_tet = new int*[ntets];
        for(int i = 0; i < ntets; i++)
        signs_tet[i] = new int[6];
        for(int i=0;i<ntets;i++){
            for(int j=0;j<6;j++){
                signs_tet[i][j]=(tetrahedrons(i,Edge_seq[j][0])-tetrahedrons(i,Edge_seq[j][1]))<0?-1:1;
            }
        }
        int FEdge_seq[3][2]={{0, 1}, {1, 2},{2 ,0}};
        signs_tri = new int*[ntris];
        for(int i = 0; i < ntris; i++)
        signs_tet[i] = new int[3];
        for(int i=0;i<ntets;i++){
            for(int j=0;j<3;j++){
                signs_tri[i][j]=(triangles(i,FEdge_seq[j][0])-triangles(i,FEdge_seq[j][1]))<0?-1:1;
            }
        }


}
void MeshEngine::getEdgeData(){
    list <int> lst;
    int id;
    int cnt=0;
    int Edge_seq[6][2]={{0, 1}, {0, 2},{0 ,3}, {1, 2}, {3 ,1}, {2, 3}};
    int FEdge_seq[3][2]={{0, 1}, {1, 2},{2 ,0}};

    for(int i=0;i<ntets;i++){
        for(int j=0;j<6;j++){
            lst.push_back(getId(tetrahedrons(i,Edge_seq[j][0]),tetrahedrons(i,Edge_seq[j][1])));
        }
    }

    lst.sort();
    lst.unique();


    tet2edge = new int*[ntets];
    for(int i = 0; i < ntets; i++)
     tet2edge[i] = new int[6];
 nedges=(int)lst.size();

 int *idarray = new int[nedges];
 copy(lst.begin(),lst.end(),idarray);
 for(int i=0;i<ntets;i++){

    for(int j=0;j<6;j++){


        id=getId(tetrahedrons(i,Edge_seq[j][0]),tetrahedrons(i,Edge_seq[j][1]));
        tet2edge[i][j]= Util::binarySearch(idarray,0,static_cast<int>(lst.size()),id);

    }
}
tri2edge.zeros(ntris,3);

for(int i=0;i<ntris;i++){
    for(int j=0;j<3;j++){
        id=getId(triangles(i,FEdge_seq[j][0]),triangles(i,FEdge_seq[j][1]));
        tri2edge.set(i,j, Util::binarySearch(idarray,0,static_cast<int>(lst.size()),id));
    }
}

}
 int MeshEngine::getId(int a, int b){
    return static_cast<int>(a<b?0.5*(a+b)*(a+b+1)+b:0.5*(a+b)*(a+b+1)+a);
}

void MeshEngine::getMeshInfo(){
    string line;
    ifstream meshfile(meshfilename+".mphtxt");
    while(getline(meshfile,line)){
        if(line.find("# number of mesh points")!= string::npos){
            vector <string> temp = Util::split(line,"#");
            nnodes=stoi(temp[0]);
            nodes.zeros(nnodes,3);
            
               // nodes = (double*)malloc(3*nnodes*sizeof(double));

            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);

            int ix=0;
            while(getline(meshfile,line)){

                if(line.empty()) break;
                temp = Util::split(line," ");
                    //cout<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;
                //cout<<stod(temp[0])<<endl;
                nodes.set(ix,0,stod(temp[0]));
                nodes.set(ix,1,stod(temp[1]));
                nodes.set(ix,2,stod(temp[2]));
                ix++;
            }
            
          
            
        }
        else if(line=="# Type #2"){
                //cout<<"triangle";
            
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);

            vector <string> temp = Util::split(line,"#");

            ntris=stoi(Util::trim(temp[0]));

            triangles.zeros(ntris,3);

            tritype.zeros(ntris,1);
            getline(meshfile,line);
            int ix=0;

            while(getline(meshfile,line)){

                if(line.empty()) break;
                temp = Util::split(line," ");
                   // cout<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;
                    //*(triangles+ix)=stoi(temp[0]);
                    //*(triangles+ntris+ix)=stoi(temp[1]);
                    //*(triangles+ntris*2+ix)=stoi(temp[2]);
                triangles.set(ix,0,stoi(temp[0]));
                triangles.set(ix,1,stoi(temp[1]));
                triangles.set(ix,2,stoi(temp[2]));
                ix++;
            }

            getline(meshfile,line);
            getline(meshfile,line);
            
            ix=0;




            int tt;

            while(getline(meshfile,line)){

                if(line.empty()) break;
                tt=stoi(Util::trim(line));
               
                if(boundDictionary[tt+1]!=0){
                    tritype[ix]=boundDictionary[tt+1];

                }
                else{
                    tritype[ix]=-1;
                
             }
                ix++;
            }
            


            
        }else if(line=="# Type #3"){

            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            getline(meshfile,line);
            vector <string> temp = Util::split(line,"#");
            ntets=stoi(temp[0]);
               // cout<<ntets;
            tetrahedrons.zeros(ntets,4);
            
            getline(meshfile,line);
            int ix=0;
            tettype= new int[ntets];

            while(getline(meshfile,line)){

                if(line.empty()) break;
                temp = Util::split(line," ");
                  // cout<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;
                tetrahedrons.set(ix,0,stoi(temp[0]));
                tetrahedrons.set(ix,1,stoi(temp[1]));
                tetrahedrons.set(ix,2,stoi(temp[2]));
                tetrahedrons.set(ix,3,stoi(temp[3]));
                ix++;
            }

      
            getline(meshfile,line);
            getline(meshfile,line);
            //getline(meshfile,line);
            
            ix=0;
            while(getline(meshfile,line)){
                
                //cout<<domDictionary[stoi(line)];

                if(Util::trim(line).empty()) {cout<<"breaking";break;}

                tettype[ix]=domDictionary[stoi(line)];

                ix++;
            }
           
            return;
                //cout<<ix;
        }
    }

}
void MeshEngine::getTagInfo(){

    pipeline pl=pipeline(24);
    string str;
    int mode;
    mode=1;
    int modecount=0;
    int nmode=0;

    ifstream meshfile(meshfilename+".html");

    while(mode!=100){
        
        while(getline(meshfile,str)){
            
            pl.push(str);
            if(str.find("Selections</span>")!= string::npos){
                nmode=1;
                modecount++;

            }

            if(str.find("Coordinate Systems</span>")!= string::npos){

                nmode=0;


            }

            

            if(str.find("Boundaries")!= string::npos){
                if(nmode==1){

                    mode=2;
                   // nmode=0;
                    modecount++;
                    break;
                }
            }

            if(str.find("Boundary")!= string::npos){
                if(nmode==1){

                       // nmode=0;
                    mode=3;
                    modecount++;
                    break;
                }
            }

            if(str.find("Domains")!= string::npos){
               if(nmode==1){

                mode=4;
                modecount++;
                break;
            }
        }

        if(str.find("Domain")!= string::npos && str.find("Domains")== string::npos){
            if(nmode==1){

                mode=5;
                modecount++;

                break;
            }
        }

        if(str.find("<span>Equations</span>")!= string::npos){

            mode=100;
            modecount++;
            break;
        }

        if(str.find("</html>")!= string::npos){
            mode=100;
            modecount++;
            break;
        }

    }
    if(mode==2){
        
        TagSet ts;
        string temp3=pl.fetch();
        temp3.replace( temp3.find("</span>"), 8," ");
        temp3.replace( temp3.find("<span>"), 6," ");
        ts.name = Util::trim(temp3);
        if(ts.name=="IP"){
            ts.tagId=2;
        }else if(ts.name=="OP"){
            ts.tagId=3;
        }else if(ts.name=="PEC"){
            ts.tagId=1;
        }else{
            boundcount++;
            ts.tagId=boundcount;
        }
        replace( str.begin(), str.end(), (char) -30, ' ');
        replace( str.begin(), str.end(), (char) -128, ' ');
        replace( str.begin(), str.end(), (char) -109, '-');
        str.replace( str.find("<span>Boundaries"), 16," ");
        str.replace( str.find("</span>"), 8," ");
        vector<string> temp= Util::split(str,",");

        for(int i=0;i<temp.size();i++)
        {  
            vector<string> temp2= Util::split(temp[i],"-");
            if(temp2.size()==2){
                for(int j=stoi(Util::trim(temp2[0]));j<=stoi(Util::trim(temp2[1]));j++)
                    {ts.tags.push_back(j);
                        boundDictionary.insert(pair <int,int> (j,ts.tagId));}
                    }else{
                       ts.tags.push_back(stoi(Util::trim(temp2[0])));
                   }
               }

               Boundaries.push_back(ts);




           }else if(mode==3){
            TagSet ts;

            str.replace( str.find("<span>Boundary"), 14," ");
            str.replace( str.find("</span>"), 8," ");

            ts.tags.push_back(stoi(Util::trim(str)));
            string temp3=pl.fetch();
            temp3.replace( temp3.find("</span>"), 8," ");
            temp3.replace( temp3.find("<span>"), 6," ");
            ts.name = Util::trim(temp3);
            if(ts.name=="IP"){
                ts.tagId=2;
            }else if(ts.name=="OP"){
                ts.tagId=3;
            }else if(ts.name=="PEC"){
                ts.tagId=1;
            }else{
                boundcount++;
                ts.tagId=boundcount;
            }
            Boundaries.push_back(ts);
            boundDictionary.insert(pair <int,int> (stoi(Util::trim(str)),ts.tagId));
        }else if(mode==4){
            TagSet ts;


            replace( str.begin(), str.end(), (char) -30, ' ');
            replace( str.begin(), str.end(), (char) -128, ' ');
            replace( str.begin(), str.end(), (char) -109, '-');
            str.replace( str.find("<span>Domains"), 13," ");
            str.replace( str.find("</span>"), 8," ");

            vector<string> temp= Util::split(str,",");

            for(int i=0;i<temp.size();i++)
            {   

                vector<string> temp2= Util::split(temp[i],"-");
                if(temp2.size()==2){
                    for(int j=stoi(Util::trim(temp2[0]));j<=stoi(Util::trim(temp2[1]));j++){

                        ts.tags.push_back(j);
                        domDictionary.insert(pair <int,int> (j,ts.tagId));
                    }

                }else{
                   ts.tags.push_back(stoi(Util::trim(temp2[0])));
               }
           }

           string temp3=pl.fetch();


           temp3.replace( temp3.find("</span>"), 8," ");
           temp3.replace( temp3.find("<span>"), 6," ");
           ts.name = Util::trim(temp3);
           if(ts.name=="AIR"){
            ts.tagId=0;
        }else{
            domcount++;;
            ts.tagId=domcount;;
        }

        Domains.push_back(ts);





    }else if(mode==5){
       TagSet ts;

       str.replace( str.find("<span>Domain"), 12," ");
       str.replace( str.find("</span>"), 8," ");


       ts.tags.push_back(stoi(Util::trim(str)));
       string temp3=pl.fetch();
       
       
       temp3.replace( temp3.find("</span>"),8," ");
      
       temp3.replace( temp3.find("<span>"), 6," ");
       ts.name = Util::trim(temp3);

       if(ts.name=="AIR"){
        ts.tagId=0;
    }else{
        domcount++;;
        ts.tagId=domcount;;
    }
    Domains.push_back(ts);
    domDictionary.insert(pair <int,int> (stoi(Util::trim(str)),ts.tagId));
}
}


meshfile.close();
}
