// dns.cpp
#include "DNS.h"

// 全局变量定义
int n_x0, n_y0;
double dx, dy, vx;
double velocity;
double l2_norm_x = 0.0, l2_norm_y = 0.0, l2_norm_p = 0.0;
double a, b;

// Mesh 类的构造函数
Mesh::Mesh(int n_y, int n_x)
    : u(n_y + 2, n_x + 2), u_star(n_y + 2, n_x + 2),u0(n_y + 2, n_x + 2), du(n_y + 2, n_x + 2),
      v(n_y + 2, n_x + 2), v_star(n_y + 2, n_x + 2),v0(n_y + 2, n_x + 2), dv(n_y + 2, n_x + 2),
      p(n_y + 2, n_x + 2), p_star(n_y + 2, n_x + 2), p_prime(n_y + 2, n_x + 2),
      u_face(n_y + 2, n_x + 1), v_face(n_y + 1, n_x + 2) {}

// 初始化所有矩阵为零
void Mesh::initializeToZero() {
    u.setZero();
    u_star.setZero();
    v.setZero();
    v_star.setZero();
    p.setZero();
    p_star.setZero();
    p_prime.setZero();
    u_face.setZero();
    v_face.setZero();
    u0.setZero();
    du.setZero();
    v0.setZero();
    dv.setZero();
}

// 显示矩阵内容
void Mesh::displayMatrix(const Eigen::MatrixXd& matrix, const std::string& name) const {
    std::cout << name << ":\n" << matrix << "\n";
}

// 显示所有矩阵
void Mesh::displayAll() const {
    displayMatrix(u, "u");
    displayMatrix(u_star, "u_star");
    displayMatrix(v, "v");
    displayMatrix(v_star, "v_star");
    displayMatrix(p, "p");
    displayMatrix(p_star, "p_star");
    displayMatrix(p_prime, "p_prime");
    displayMatrix(u_face, "u_face");
    displayMatrix(v_face, "v_face");
}

// Equation 类的构造函数
Equation::Equation(int n_y_, int n_x_)
    : A_p(n_y_ + 2, n_x_ + 2), A_e(n_y_ + 2, n_x_ + 2), A_w(n_y_ + 2, n_x_ + 2),
      A_n(n_y_ + 2, n_x_ + 2), A_s(n_y_ + 2, n_x_ + 2), source(n_y_ * n_x_),
      A(n_y_ * n_x_, n_y_ * n_x_), n_x(n_x_), n_y(n_y_) {}

// 初始化矩阵和源向量为零
void Equation::initializeToZero() {
    A_p.setZero();
    A_e.setZero();
    A_w.setZero();
    A_n.setZero();
    A_s.setZero();
    source.setZero();
    A.setZero();
}

// 使用输入的系数矩阵构建稀疏矩阵 A
void Equation::build_matrix() {
            typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        int n;

        // 顶盖
        int i = 1;
        for (int j = 2; j < n_x; j++) {
            n = (j - 1);
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        //内部
        for (i = 2; i < n_y; i++) {
            for (int j = 2; j < n_x; j++) {
                n = (i - 1) * n_x + (j - 1);
                tripletList.emplace_back(n, n, A_p(i, j));
                tripletList.emplace_back(n, n - 1, -A_w(i, j));
                tripletList.emplace_back(n, n + 1, -A_e(i, j));
                tripletList.emplace_back(n, n - n_x, -A_n(i, j));
                tripletList.emplace_back(n, n + n_x, -A_s(i, j));
            }
        }

        int j = 1;
        for (i = 2; i < n_y; i++) {
            n = (i - 1) * n_x;
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        j = n_x;
        for (i = 2; i < n_y; i++) {
            n = i * n_x - 1;
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        i = n_y;
        for (j = 2; j < n_x; j++) {
            n = (n_y - 1) * n_x + (j - 1);
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
        }

        
        n = 0;
        i = 1;
        j = 1;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n + 1, -A_e(i, j));
        tripletList.emplace_back(n, n + n_x, -A_s(i, j));

        
        i = 1;
        j = n_x;
        n = n_x - 1;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n - 1, -A_w(i, j));
        tripletList.emplace_back(n, n + n_x, -A_s(i, j));

        
        n = (n_y - 1) * n_x;
        j = 1;
        i = n_y;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n + 1, -A_e(i, j));
        tripletList.emplace_back(n, n - n_x, -A_n(i, j));

        
        n = (n_x * n_y) - 1;
        i = n_y;
        j = n_x;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n - 1, -A_w(i, j));
        tripletList.emplace_back(n, n - n_x, -A_n(i, j));

        
        A.setFromTriplets(tripletList.begin(), tripletList.end());
    }
void solve(SparseMatrix <double> &A,VectorXd &b,MatrixXd &phi,double &l2_norm,double &epsilon,int n_x,int n_y)
{
    
    VectorXd x(n_x*n_y);

    int n=0;
    for(int i=1;i<n_y+1;i++)
    {
        for (int j=1;j<n_x+1;j++)
        {
         
         x[n]=phi(i,j);
         n+=1;

        }
    }
    l2_norm=(A*x-b).norm() ;

    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(A);
    solver.setTolerance(epsilon);
    x = solver.solve(b);
    n=0;
    for(int i=1;i<n_y+1;i++)
    {
        for (int j=1;j<n_x+1;j++)
        {
         
         phi(i,j)=x[n];
         n+=1;

        }
    }



}


void face_velocity(Mesh &mesh,Equation &equ_u)
{
    MatrixXd &u_face= mesh.u_face;
    MatrixXd &v_face= mesh.v_face;
    MatrixXd &u= mesh.u;
    MatrixXd &v= mesh.v;
    MatrixXd &p= mesh.p;
    MatrixXd &p_star= mesh.p_star;
    MatrixXd &p_prime= mesh.p_prime;
    MatrixXd &u_star= mesh.u_star;
    MatrixXd &v_star= mesh.v_star;
    MatrixXd &A_p=equ_u.A_p;
    MatrixXd &A_e=equ_u.A_e;
    MatrixXd &A_w=equ_u.A_w;
    MatrixXd &A_n=equ_u.A_n;
    MatrixXd &A_s=equ_u.A_s;
    SparseMatrix<double> &A=equ_u.A;
    int n_y=equ_u.n_y;
    int n_x=equ_u.n_x;
    //x面上速度
    int i,j;
    double alpha_uv=10e-2;
    for(i=1;i<n_y+1;i++)
    {
      for (j=1;j<n_x;j++)
      {
        u_face(i,j)=0.5*(u(i,j) + u(i,j + 1)) + 0.25*alpha_uv*(p(i,j + 1) - p(i,j - 1))*dy/A_p(i,j) + 0.25*alpha_uv*(p(i,j + 2) - p(i,j))*dy/A_p(i,j + 1)- 0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(p(i,j + 1) - p(i,j))*dy;
          
      }  
    }



    //y面上熟读
    for(i=2;i<n_y+1;i++)
    {
        for (j=1;j<n_x+1;j++)
        {

            v_face(i-1,j)=0.5*(v(i,j) + v(i - 1,j)) + 0.25*alpha_uv*(p(i - 1,j) - p(i + 1,j))*dy/A_p(i,j) + 0.25*alpha_uv*(p(i - 2,j) - p(i,j))*dy/A_p(i - 1,j)- 0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(p(i - 1,j) - p(i,j))*dy;
        }
    }

}


void pressure_function(Mesh &mesh,Equation &equ_p,Equation &equ_u)
{
    MatrixXd &u_face= mesh.u_face;
    MatrixXd &v_face= mesh.v_face;
    MatrixXd &u= mesh.u;
    MatrixXd &v= mesh.v;
    MatrixXd &p= mesh.p;
    MatrixXd &p_star= mesh.p_star;
    MatrixXd &p_prime= mesh.p_prime;
    MatrixXd &u_star= mesh.u_star;
    MatrixXd &v_star= mesh.v_star;
    MatrixXd &Ap_p=equ_p.A_p;
    MatrixXd &Ap_e=equ_p.A_e;
    MatrixXd &Ap_w=equ_p.A_w;
    MatrixXd &Ap_n=equ_p.A_n;
    MatrixXd &Ap_s=equ_p.A_s;
    MatrixXd &A_p=equ_u.A_p;
    Eigen::VectorXd &source_p=equ_p.source;
    int n_x=equ_u.n_x;
    int n_y=equ_u.n_y;
    int n,i,j;
    double alpha_uv=10e-2;
    //i内部面
    for(i=2;i<n_y;i++)
    {
        for(j=2;j<n_x;j++)
        {
            n=(i-1)*n_x + (j-1);

            Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
            Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
            Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
            Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
            Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

            source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;

        } 

    }
    
    //顶面
    i=1;
    for(j=2;j<n_x;j++)
    { 
        n=(j-1);
    
        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }
    
    //左壁面
    j=1;
    for(i=2;i<n_y;i++)
    { 
        n=(i-1)*n_x ;

        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }

    //有壁面
    j=n_x;
    for(i=2;i<n_y;i++)
    {
        n=i*n_x -1;    

        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }
    //bottom
    i=n_y;
    for(j=2;j<n_x;j++)
    { 
        n=(n_y-1)*n_x + (j-1);

        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }

    //top left corner
    i=1;
    j=1;
    n=0;

    Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
    Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
    //top right corner
    i=1;
    j=n_x;
    n=n_x-1;

    Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
    Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;

    i=n_y;
    j=1;
    n=(n_y-1)*n_x ;

    Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
    Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
 
    i=n_y;
    j=n_x;
    n=n_x*n_y-1;

    Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
    Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
}

//修正压力
void correct_pressure(Mesh &mesh,Equation &equ_u)
{
    MatrixXd &p= mesh.p;
    MatrixXd &p_star= mesh.p_star;
    MatrixXd &p_prime= mesh.p_prime;
    int n_x=equ_u.n_x;
    int n_y=equ_u.n_y;
    p_prime.block(0,1,1,n_x)=p_prime.block(1,1,1,n_x);
    
   
    p_prime.block(1,0,n_y,1)=p_prime.block(1,1,n_y,1);
    

    p_prime.block(1,n_x+1,n_y,1)=p_prime.block(1,n_x,n_y,1);
    
    p_prime.block(n_y+1,1,1,n_x)=p_prime.block(n_y,1,1,n_x);
    

    p_prime(0,0)=(p_prime(1,1)+p_prime(0,1)+p_prime(1,0))/3;

    
    p_prime(0,n_x+1)=(p_prime(0,n_x)+p_prime(1,n_x)+p_prime(1,n_x+1))/3;

    p_prime(n_y+1,0)=(p_prime(n_y,0)+p_prime(n_y,1)+p_prime(n_y+1,1))/3;

 
    p_prime(n_y+1,n_x+1)=(p_prime(n_y,n_x+1)+p_prime(n_y+1,n_x)+p_prime(n_y,n_x))/3;

    MatrixXd p_ref=MatrixXd::Constant(n_y+2,n_x+2,p_prime(0,0));
    
    double alpha_p=0.1;
    p_star=p+alpha_p*(p_prime);



}


void correct_velocity(Mesh &mesh,Equation &equ_u)

{
    MatrixXd &u= mesh.u;
    MatrixXd &v= mesh.v;
    MatrixXd &u_face= mesh.u_face;
    MatrixXd &v_face= mesh.v_face;
    MatrixXd &p= mesh.p;
    MatrixXd &p_star= mesh.p_star;
    MatrixXd &p_prime= mesh.p_prime;
    MatrixXd &u_star= mesh.u_star;
    MatrixXd &v_star= mesh.v_star;
    MatrixXd &A_p=equ_u.A_p;
    MatrixXd &A_e=equ_u.A_e;
    MatrixXd &A_w=equ_u.A_w;
    MatrixXd &A_n=equ_u.A_n;
    MatrixXd &A_s=equ_u.A_s;
    SparseMatrix<double> &A=equ_u.A;
    int i,j;
    int n_x=equ_u.n_x;
    int n_y=equ_u.n_y;
    double alpha_uv=10e-2;
    //u 速度修正
    //内部面
    for(i=1;i<n_y+1;i++) 
    {
        for(j=2;j<n_x;j++) 
        {
        
            u_star(i,j)= u(i,j) + 0.5*alpha_uv*(p_prime(i,j-1)-p_prime(i,j+1))*dy/A_p(i,j);
        }
    }


    //左壁面
    j=1;
    for(i=1;i<n_y+1;i++)
    {
        u_star(i,j)=u(i,j) + 0.5*alpha_uv*(p_prime(i,j) - p_prime(i,j+1))*dy/A_p(i,j);
    }

    //右壁面
    j=n_x;
    for(i=1;i<n_y+1;i++)
    {
        u_star(i,j)=u(i,j) + 0.5*alpha_uv*(p_prime(i,j-1) - p_prime(i,j))*dy/A_p(i,j);
    }

    //v速度修正
    //内部面修正
    for(i=2;i<n_y;i++) 
    {
        for(j=1;j<n_x+1;j++) 
        {
            v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i+1,j)-p_prime(i-1,j))*dx/A_p(i,j);
        }
    }


    //顶面
    i=1;
    for(j=1;j<n_x+1;j++)
    { 
        v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i + 1,j) - p_prime(i,j))*dx/A_p(i,j);
    }

    //底面
    i=n_y;
    for(j=1;j<n_x+1;j++)
    {
        v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i,j) - p_prime(i - 1,j))*dx/A_p(i,j);
    }
    for(i=1;i<n_y+1;i++)
    { 
        for(j=1;j<n_x;j++)
        {
            u_face(i,j)=u_face(i,j)+ 0.5*alpha_uv*(1/A_p(i,j)+1/A_p(i,j+1))*(p_prime(i,j)-p_prime(i,j+1))*dy;
        }
    }

    for(i=2;i<n_y+1;i++)
    {
        for(j=1;j<n_x+1;j++) 
        {
            v_face(i-1,j)=v_face(i-1,j) +  0.5*alpha_uv*(1/A_p(i,j)+1/A_p(i-1,j))*(p_prime(i,j)-p_prime(i-1,j))*dx;
        }
    }


}


void post_processing(Mesh &mseh,int n_x,int n_y,double a)
{   
    VectorXd x(n_x+2),y(n_y+2);
    x << 0,VectorXd::LinSpaced(n_x,dx/2.0,a-dx/2.0),a;
    y << 0,VectorXd::LinSpaced(n_y,dy/2.0,a-dy/2.0),a
    ;

    //保存计算结果
     std::ofstream outFile;
     outFile.open("u.dat");
     outFile << mseh.u_star;
     outFile.close();

     outFile.open("v.dat");
     outFile << mseh.v_star;
     outFile.close();

     outFile.open("p.dat");
     outFile << mseh.p_star;
     outFile.close();

     outFile.open("x.dat");
     outFile << x;
     outFile.close();

     outFile.open("y.dat");
     outFile << y;
     outFile.close();




}