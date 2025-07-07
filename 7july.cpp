#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include<string>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/SparseLU>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include<set>
#include<tuple>
using namespace std;
using namespace Eigen;
typedef Triplet<double> Tri;
vector<vector<double>> readCSV(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    string line;
    getline(file,line);
    while (getline(file,line)) {
        stringstream ss(line);
        string val;
        vector<double> row;
        while (getline(ss,val,',')) {
            row.push_back(stod(val));
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    return data;
}
struct point{
    double x,y,z;
    int voxel; 
    point(double i, double j, double k) {
        x=i, y=j, z=k;
        voxel=0;
    }
    bool operator<(const point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};
int find_voxel(double x, double y, double z, double x_min, double x_max, double y_min, double y_max,double z_min, double z_max, int voxels_inrow, int voxels_incolumn,double radius) {
    return int((z-z_min)/radius)*voxels_incolumn*voxels_inrow+int((y-y_min)/radius)*voxels_inrow+int((x-x_min)/radius);
}
int provide_znumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return int(voxel/(voxels_inrow*voxels_incolumn));
}
int provide_ynumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return (voxel-provide_znumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow*voxels_incolumn)/voxels_inrow;
}
int provide_xnumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return (voxel-provide_znumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow*voxels_incolumn-provide_ynumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow);
}
double gaussian_weight_function(point Ni, point p0, double radius) {
    double distance=(pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2)+pow(Ni.z-p0.z,2))/(radius*radius);
    if (distance<=1) {
        return exp(-6.25*distance);
    }
    else {
        return 0.0;
    }
}
int main() {
    //taking data from csv file
    vector<vector<double>> dataofvessel=readCSV("vessel_points.csv");
    vector<vector<double>> dataofbranch=readCSV("branch_points.csv");
    double x_max=INT_MIN, x_min=INT_MAX, y_max=INT_MIN, y_min=INT_MAX, z_min=INT_MAX, z_max=INT_MIN;
    for (int i=0; i<dataofvessel.size(); i++) {
        x_max=max(x_max,dataofvessel[i][0]);
        x_min=min(x_min,dataofvessel[i][0]);
        y_max=max(y_max,dataofvessel[i][1]);
        y_min=min(y_min,dataofvessel[i][1]);
        z_min=min(z_min,dataofvessel[i][2]);
        z_max=max(z_max,dataofvessel[i][2]);
    }
    for (int i=0; i<dataofbranch.size(); i++) {
        x_max=max(x_max,dataofbranch[i][0]);
        x_min=min(x_min,dataofbranch[i][0]);
        y_max=max(y_max,dataofbranch[i][1]);
        y_min=min(y_min,dataofbranch[i][1]);
        z_min=min(z_min,dataofbranch[i][2]);
        z_max=max(z_max,dataofbranch[i][2]);
    }
    //parameters of blood
    double mew=3.5*pow(10,-6), rhob=1058;
    double radius=0.0005*3.0;
    int voxels_inrow=int((x_max-x_min)/radius)+1;
    int voxels_incolumn=int((y_max-y_min)/radius)+1;
    int voxels_inz=int((z_max-z_min)/radius)+1;
    int max_voxel=voxels_inrow*voxels_incolumn*voxels_inz;
    vector<point> domain;
    map<point,int> type;
    map<point,string> vesselorbranch;
    for (int i=0; i<dataofvessel.size(); i++) {
        point p0=point(dataofvessel[i][0],dataofvessel[i][1],dataofvessel[i][2]);
        p0.voxel=find_voxel(p0.x,p0.y,p0.z,x_min,x_max,y_min,y_max,z_min,z_max,voxels_inrow,voxels_incolumn,radius);
        type[p0]=dataofvessel[i][3];
        vesselorbranch[p0]="vessel";
        domain.push_back(p0);
    }
    for (int i=0; i<dataofbranch.size(); i++) {
        point p0=point(dataofbranch[i][0],dataofbranch[i][1],dataofbranch[i][2]);
        p0.voxel=find_voxel(p0.x,p0.y,p0.z,x_min,x_max,y_min,y_max,z_min,z_max,voxels_inrow,voxels_incolumn,radius);
        type[p0]=dataofbranch[i][3];
        vesselorbranch[p0]="branch";
        domain.push_back(p0);
    }
    auto point_key = [](const point& p) {
    // Round to 6 decimal places
    double rx = round(p.x * 1e6) / 1e6;
    double ry = round(p.y * 1e6) / 1e6;
    double rz = round(p.z * 1e6) / 1e6;
    return std::make_tuple(rx, ry, rz);
};

std::set<std::tuple<double, double, double>> seen;
std::vector<point> updated_domain;

for (const point& p : domain) {
    auto key = point_key(p);
    if (seen.find(key) == seen.end()) {
        seen.insert(key);
        updated_domain.push_back(p);
    }
}

std::cout << "Original domain size: " << domain.size() << std::endl;
std::cout << "Unique domain size:   " << updated_domain.size() << std::endl;
std::cout << "Duplicates removed:   " << domain.size() - updated_domain.size() << std::endl;

// Replace old domain with the unique one
domain = updated_domain;
    vector<vector<point>> points_insidevoxel(max_voxel+1);
    for (int i=0; i<domain.size(); i++) {
        point p0=domain[i];
        points_insidevoxel[p0.voxel].push_back(p0);
    }
    vector<vector<int>> neighbours_ofvoxel(max_voxel+1);
    for(int i=0; i<=max_voxel; i++) {
        int find_x=provide_xnumber(i,voxels_inrow,voxels_incolumn);
        int find_y=provide_ynumber(i,voxels_inrow,voxels_incolumn);
        int find_z=provide_znumber(i,voxels_inrow,voxels_incolumn);
        for (int diffx=-1; diffx<=1; diffx++) {
            for (int diffy=-1; diffy<=1; diffy++) {
                for (int diffz=-1; diffz<=1; diffz++) {
                    if (find_x+diffx>=0 && find_x+diffx<=voxels_inrow-1 && find_y+diffy>=0 && find_y+diffy<=voxels_incolumn-1 && diffz+find_z>=0 && find_z+diffz<=voxels_inz-1) {
                            neighbours_ofvoxel[i].push_back((find_z+diffz)*voxels_incolumn*voxels_inrow+(find_y+diffy)*voxels_inrow+(find_x+diffx));
                }
                }
            }
        }
    }
    cout<<max_voxel<<endl;
    map<point,vector<point>> neighbours_ofpoint;
    for (int i=0; i<domain.size(); i++) {
        point p0=domain[i];
        int voxel_number=p0.voxel;
        for (int i=0; i<neighbours_ofvoxel[voxel_number].size(); i++) {
            int neighbour_voxel=neighbours_ofvoxel[voxel_number][i];
            for (int j=0; j<points_insidevoxel[neighbour_voxel].size(); j++) {
                point Ni=points_insidevoxel[neighbour_voxel][j];
                double distance=pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2)+pow(Ni.z-p0.z,2);
                if (distance>0 && distance<=pow(radius,2)) {
                    neighbours_ofpoint[p0].push_back(Ni);
                }
                if (neighbours_ofpoint[p0].size()>=100) {
                    break;
                }
            }
            if (neighbours_ofpoint[p0].size()>=100) {
                    break;
                }
        }
    }
    int max_neighbours=INT_MIN, min_neighbours=INT_MAX;
    for (int i=0; i<domain.size(); i++) {
        point p0=domain[i];
        int k=neighbours_ofpoint[p0].size();
        max_neighbours=max(max_neighbours,k);
        min_neighbours=min(min_neighbours,k);
}
    cout<<"maximum neighbours: "<<max_neighbours<<" minimum neighbours: "<<min_neighbours<<endl;
    //declaring previous velocities
    map<point,double> prev_vx,prev_vy,prev_vz;
    for (int p=0; p<domain.size(); p++) {
        point p0=domain[p];
        if (type[p0]==0 || type[p0]==1) {
            prev_vx[p0]=0.0, prev_vy[p0]=0.0, prev_vz[p0]=0.0;
        }
        else if (type[p0]==2) {
            prev_vx[p0]=0.0, prev_vy[p0]=0.0, prev_vz[p0]=0.3;
        }
        else if (type[p0]==3) {
            prev_vx[p0]=0.15, prev_vy[p0]=0.0, prev_vz[p0]=0.0;
        }
        else if (type[p0]==4) {
            prev_vx[p0]=0.0, prev_vy[p0]=0.0, prev_vz[p0]=0.2625;
        }
    }
    //starting of a time loop
    double T,dt;
    cout<<"Enter the maximum time: ";
    cin>>T;
    cout<<"Enter the time step: ";
    cin>>dt;
    int total_steps=int(T/dt);
    for (int t=1; t<=total_steps; t++) {
        cout<<"time step- "<<t<<" started"<<endl;
        //calculating squarevx, squarevy, squarevz, vxvy, vyvz, vzvx
        map<point,double> square_vx, square_vy, square_vz, vxvy, vyvz, vzvx;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            square_vx[p0]=pow(prev_vx[p0],2.0);
            square_vy[p0]=pow(prev_vy[p0],2.0);
            square_vz[p0]=pow(prev_vz[p0],2.0);
            vxvy[p0]=prev_vx[p0]*prev_vy[p0];
            vyvz[p0]=prev_vy[p0]*prev_vz[p0];
            vzvx[p0]=prev_vz[p0]*prev_vx[p0];
        }
        //calculating laplacian of vx, vy, vz
        map<point,double> laplacian_vx, laplacian_vy, laplacian_vz;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,10);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b1(totalneighbours), b2(totalneighbours), b3(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                M(i,8)=dy*dz, M(i,9)=dz*dx;
                b1(i)=prev_vx[Ni];
                b2(i)=prev_vy[Ni];
                b3(i)=prev_vz[Ni];
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
            }
            VectorXd L(10);
            L<<0,0,0,0,1,1,1,0,0,0;
            VectorXd rhs1=M.transpose()*W*b1;
            MatrixXd MTWM=M.transpose()*W*M;
            VectorXd a1=MTWM.ldlt().solve(rhs1);
            laplacian_vx[p0]=L.transpose()*a1;
            VectorXd rhs2=M.transpose()*W*b2;
            VectorXd a2=MTWM.ldlt().solve(rhs2);
            laplacian_vy[p0]=L.transpose()*a2;
            VectorXd rhs3=M.transpose()* W*b3;
            VectorXd a3=MTWM.ldlt().solve(rhs3);
            laplacian_vz[p0]=L.transpose()*a3;
        }
        map<point,double> vxsquare_x, vxvy_y, vxvz_z, vxvy_x, vysquare_y, vyvz_z, vxvz_x, vyvz_y, vzsquare_z;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,10);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b1(totalneighbours), b2(totalneighbours), b3(totalneighbours),
            b4(totalneighbours), b5(totalneighbours), b6(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                M(i,8)=dy*dz, M(i,9)=dz*dx;
                b1(i)=square_vx[Ni];
                b2(i)=square_vy[Ni];
                b3(i)=square_vz[Ni];
                b4(i)=vxvy[Ni];
                b5(i)=vyvz[Ni];
                b6(i)=vzvx[Ni];
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
            }
            VectorXd L1(10), L2(10), L3(10);
            L1<<0,1,0,0,0,0,0,0,0,0;
            L2<<0,0,1,0,0,0,0,0,0,0;
            L3<<0,0,0,1,0,0,0,0,0,0;
            VectorXd rhs1=M.transpose()*W*b1;
            VectorXd rhs2=M.transpose()*W*b2;
            VectorXd rhs3=M.transpose()*W*b3;
            VectorXd rhs4=M.transpose()*W*b4;
            VectorXd rhs5=M.transpose()*W*b5;
            VectorXd rhs6=M.transpose()*W*b6;
            MatrixXd MTWM=M.transpose()*W*M;
            VectorXd a1=MTWM.ldlt().solve(rhs1);
            VectorXd a2=MTWM.ldlt().solve(rhs2);
            VectorXd a3=MTWM.ldlt().solve(rhs3);
            VectorXd a4=MTWM.ldlt().solve(rhs4);
            VectorXd a5=MTWM.ldlt().solve(rhs5);
            VectorXd a6=MTWM.ldlt().solve(rhs6);
            vxsquare_x[p0]=L1.transpose()*a1;
            vxvy_y[p0]=L2.transpose()*a4;
            vxvz_z[p0]=L3.transpose()*a6;
            vxvy_x[p0]=L1.transpose()*a4;
            vysquare_y[p0]=L2.transpose()*a2;
            vyvz_z[p0]=L3.transpose()*a5;
            vxvz_x[p0]=L1.transpose()*a6;
            vyvz_y[p0]=L2.transpose()*a5;
            vzsquare_z[p0]=L3.transpose()*a3;
        }
        //calculating v*
        map<point,double> vstar_x, vstar_y, vstar_z;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            vstar_x[p0]=prev_vx[p0]+(dt*(mew*laplacian_vx[p0]-vxsquare_x[p0]-vxvy_y[p0]-vxvz_z[p0]));
            vstar_y[p0]=prev_vy[p0]+(dt*(mew*laplacian_vy[p0]-vxvy_x[p0]-vysquare_y[p0]-vyvz_z[p0]));
            vstar_z[p0]=prev_vz[p0]+(dt*(mew*laplacian_vz[p0]-vxvz_x[p0]-vyvz_y[p0]-vzsquare_z[p0]));
        }
        cout<<"we have successfully calculated v*"<<endl;
        //calculating divergence of v*
        map<point,double> vstarx_x, vstary_y, vstarz_z;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,10);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b1(totalneighbours), b2(totalneighbours), b3(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                M(i,8)=dy*dz, M(i,9)=dz*dx;
                b1(i)=vstar_x[Ni];
                b2(i)=vstar_y[Ni];
                b3(i)=vstar_z[Ni];
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
            }
            VectorXd L1(10), L2(10), L3(10);
            L1<<0,1,0,0,0,0,0,0,0,0;
            L2<<0,0,1,0,0,0,0,0,0,0;
            L3<<0,0,0,1,0,0,0,0,0,0;
            VectorXd rhs1=M.transpose()*W*b1;
            MatrixXd MTWM=M.transpose()*W*M;
            VectorXd a1=MTWM.ldlt().solve(rhs1);
            vstarx_x[p0]=L1.transpose()*a1;
            VectorXd rhs2=M.transpose()*W*b2;
            VectorXd a2=MTWM.ldlt().solve(rhs2);
            vstary_y[p0]=L2.transpose()*a2;
            VectorXd rhs3=M.transpose()* W*b3;
            VectorXd a3=MTWM.ldlt().solve(rhs3);
            vstarz_z[p0]=L3.transpose()*a3;
        }
        cout<<"we have successfully calculated divergence of v*"<<endl;
        //calculating laplacian of pressure
        map<point,double> laplacian_p;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            laplacian_p[p0]=(-1*rhob*(vstarx_x[p0]+vstary_y[p0]+vstarz_z[p0]))/(dt);
        }
        cout<<"we have successfully calculated laplacian of pressure"<<endl;
        ofstream fout1("non_zero_entries.csv");
        fout1<<"X"<<","<<"Y"<<","<<"Z"<<","<<"non_zero_elements"<<"\n";
        //calculating pressure
        map<point,int> identity;
        int num=0;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            identity[p0]=num;
            num+=1;
        }
        map<point,double> pressure;
        int total_points=identity.size();
        SparseMatrix<double> letsSolve(total_points,total_points);
        VectorXd rhs(total_points);
        vector<Tri> coefficients;
        cout<<"size of domain: "<<domain.size()<<endl;
        cout<<"size of identity: "<<identity.size()<<endl;
        for (int p=0; p<domain.size(); p++) {
            int non_zero=0;
            double radius_vessel=0.002, radius_branch=0.001;
            point p0=domain[p];
            int tp=type[p0];
            string vob=vesselorbranch[p0];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            if (vob=="vessel") {
                if (tp==1) {
                    MatrixXd M(totalneighbours+2,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=p0.x/radius_vessel,M(totalneighbours+1,2)=p0.y/radius_vessel,M(totalneighbours+1,3)=0 ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0,
                    M(totalneighbours+1,6)=0,M(totalneighbours+1,7)=0,M(totalneighbours+1,8)=0,M(totalneighbours+1,9)=0;
                    W(totalneighbours+1,totalneighbours+1)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);
                }
                else if (tp==2 || tp==4) {
                    MatrixXd M(totalneighbours+2,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=0,M(totalneighbours+1,2)=0,M(totalneighbours+1,3)=1 ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0,
                    M(totalneighbours+1,6)=0,M(totalneighbours+1,7)=0,M(totalneighbours+1,8)=0,M(totalneighbours+1,9)=0;
                    W(totalneighbours+1,totalneighbours+1)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);
                }
                else if (tp==0) {
                    MatrixXd M(totalneighbours+1,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+1,totalneighbours+1);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);                      
                }
            }
            else if (vob=="branch")  {
                if (tp==1) {
                    MatrixXd M(totalneighbours+2,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=0,M(totalneighbours+1,2)=p0.y/radius_branch,M(totalneighbours+1,3)=(p0.z-0.025)/radius_branch,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0,
                    M(totalneighbours+1,6)=0,M(totalneighbours+1,7)=0,M(totalneighbours+1,8)=0,M(totalneighbours+1,9)=0;
                    W(totalneighbours+1,totalneighbours+1)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);
                }
                else if (tp==3) {
                    MatrixXd M(totalneighbours+2,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=1,M(totalneighbours+1,2)=0,M(totalneighbours+1,3)=0 ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0,
                    M(totalneighbours+1,6)=0,M(totalneighbours+1,7)=0,M(totalneighbours+1,8)=0,M(totalneighbours+1,9)=0;
                    W(totalneighbours+1,totalneighbours+1)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);
                }
                else if (tp==0) {
                    MatrixXd M(totalneighbours+1,10);
                    MatrixXd W =MatrixXd :: Zero(totalneighbours+1,totalneighbours+1);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                        M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                        M(i,8)=dy*dz, M(i,9)=dz*dx;
                        W(i,i)=gaussian_weight_function(Ni,p0,radius);
                    }
                    M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=0 ,M(totalneighbours,4)=1,M(totalneighbours,5)=1,
                    M(totalneighbours,6)=1,M(totalneighbours,7)=0,M(totalneighbours,8)=0,M(totalneighbours,9)=0;
                    W(totalneighbours,totalneighbours)=1;
                    MatrixXd MTWM=M.transpose()*W*M;
                    MatrixXd MTW=M.transpose()*W;  
                    MatrixXd A=MTWM.ldlt().solve(MTW);
                    for (int i=0; i<totalneighbours; i++) {
                        point Ni=neighbours[i];
                        coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                        if (A(0,i)!=0) {
                            non_zero+=1;
                        }
                    }
                    coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                    rhs(identity[p0])=-(A(0,totalneighbours)*laplacian_p[p0]);
                }
            }
            fout1<<p0.x<<","<<p0.y<<","<<p0.z<<","<<non_zero<<"\n";
        }
        fout1.close();
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        letsSolve.setFromTriplets(coefficients.begin(), coefficients.end());
        solver.compute(letsSolve);
        VectorXd answer = solver.solve(rhs);
        answer=answer.array()-answer.mean();
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            pressure[p0]=answer(identity[p0]);
        }
        cout<<"pressure calculated successfully"<<endl;
        //calculating gradient of p
        map<point,double> grad_px, grad_py, grad_pz;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,10);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double dx=Ni.x-p0.x, dy=Ni.y-p0.y, dz=Ni.z-p0.z;
                M(i,0)=1, M(i,1)=dx, M(i,2)=dy, M(i,3)=dz, M(i,4)=(dx*dx)/2, M(i,5)=(dy*dy)/2, M(i,6)=(dz*dz)/2, M(i,7)=(dx*dy),
                M(i,8)=dy*dz, M(i,9)=dz*dx;
                b(i)=pressure[Ni];
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
            }
            VectorXd L1(10), L2(10), L3(10);
            L1<<0,1,0,0,0,0,0,0,0,0;
            L2<<0,0,1,0,0,0,0,0,0,0;
            L3<<0,0,0,1,0,0,0,0,0,0;
            VectorXd rhs = M.transpose() * W * b;
            MatrixXd MTWM = M.transpose() * W * M;
            VectorXd a = MTWM.ldlt().solve(rhs);
            grad_px[p0] = L1.transpose() * a;
            grad_py[p0] = L2.transpose() * a;
            grad_pz[p0] = L3.transpose() * a;
        }
        cout<<"we have successfully calculated gradient of p"<<endl;
        //calculating new velocities
        map<point,double> next_vx,next_vy,next_vz;
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            if (type[p0]==1) {
                next_vx[p0]=0.0, next_vy[p0]=0.0, next_vz[p0]=0.0;
            }
            else if (type[p0]==2) {
                next_vx[p0]=0.0, next_vy[p0]=0.0, next_vz[p0]=0.3;
            }
            else if (type[p0]==3) {
                next_vx[p0]=0.15, next_vy[p0]=0.0, next_vz[p0]=0.0;
            }
            else if (type[p0]==4) {
                next_vx[p0]=0.0, next_vy[p0]=0.0, next_vz[p0]=0.2625;
            }
            else {
                next_vx[p0]=(-1*dt*grad_px[p0])/rhob+vstar_x[p0];
                next_vy[p0]=(-1*dt*grad_py[p0])/rhob+vstar_y[p0];
                next_vz[p0]=(-1*dt*grad_pz[p0])/rhob+vstar_z[p0];
            }
        }
        prev_vx=next_vx, prev_vy=next_vy, prev_vz=next_vz;
        double max_velocity=INT_MIN, min_velocity=INT_MAX;
        cout<<"new velocity calculated"<<endl;
        ofstream fout("Velocity"+to_string(t)+".csv");
        fout<<"X"<<","<<"Y"<<","<<"Z"<<","<<"Velocity"<<"\n";
        for (int p=0; p<domain.size(); p++) {
            point p0=domain[p];
            double vel=pow((pow(prev_vx[p0],2.0)+pow(prev_vy[p0],2.0)+pow(prev_vz[p0],2.0)),0.5);
            max_velocity=max(max_velocity,vel);
            min_velocity=min(min_velocity,vel);
            fout<<p0.x<<","<<p0.y<<","<<p0.z<<","<<vel<<"\n";
        }
        fout.close();
        cout<<"minimum velocity: "<<min_velocity<<" maximum velocity: "<<max_velocity<<endl;
        cout<<"time step- "<<t<<" ended"<<endl;
    }
    return 0;
}