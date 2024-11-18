#include <iostream>
#include  <eigen3/Eigen/Sparse>
#include<algorithm> 
#include<fstream>

using namespace Eigen;
    //声明并初始变量
    int n_x0,n_y0;
    double dx,dy,vx;
    double velocity;

    double l2_norm_x=0.0;
    double l2_norm_y=0.0;
    double l2_norm_p=0.0;
    
    double a ;
    double b ;
class Mesh {
public:
    // 成员变量
    Eigen::MatrixXd u, u_star;
    Eigen::MatrixXd v, v_star;
    Eigen::MatrixXd p, p_star, p_prime;
    Eigen::MatrixXd u_face, v_face;

    // 默认构造函数
    Mesh() = default;

    // 参数化构造函数
    Mesh(int n_y, int n_x)
        : u(n_y + 2, n_x + 2), u_star(n_y + 2, n_x + 2),
          v(n_y + 2, n_x + 2), v_star(n_y + 2, n_x + 2),
          p(n_y + 2, n_x + 2), p_star(n_y + 2, n_x + 2), p_prime(n_y + 2, n_x + 2),
          u_face(n_y + 2, n_x + 1), v_face(n_y + 1, n_x + 2) {}

    // 初始化所有矩阵为零
    void initializeToZero() {
        u.setZero();
        u_star.setZero();
        v.setZero();
        v_star.setZero();
        p.setZero();
        p_star.setZero();
        p_prime.setZero();
        u_face.setZero();
        v_face.setZero();
    }

    // 显示矩阵内容
    void displayMatrix(const Eigen::MatrixXd& matrix, const std::string& name) const {
        std::cout << name << ":\n" << matrix << "\n";
    }

    // 显示所有矩阵
    void displayAll() const {
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
};
class Equation {
public:
    // Member variables
    Eigen::MatrixXd A_p, A_e, A_w, A_n, A_s;
    Eigen::VectorXd source;
    Eigen::SparseMatrix<double> A; // Sparse matrix for coefficients
    int n_x, n_y;                 // Grid dimensions

    // Constructor
    Equation(int n_y_, int n_x_)
        : A_p(n_y_+2, n_x_+2), A_e(n_y_+2, n_x_+2), A_w(n_y_+2, n_x_+2), 
          A_n(n_y_+2, n_x_+2), A_s(n_y_+2, n_x_+2), source(n_y_ * n_x_), 
          A((n_y_) *( n_x_), (n_y_) *( n_x_)), n_x(n_x_), n_y(n_y_) {}

    // Initialize matrices and source vector to zero
    void initializeToZero() {
        A_p.setZero();
        A_e.setZero();
        A_w.setZero();
        A_n.setZero();
        A_s.setZero();
        source.setZero();
        A.setZero();
    }

    // Method to build the sparse matrix A using input coefficient matrices
    void build_matrix() {
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        int n;

        // Top boundary
        int i = 1;
        for (int j = 2; j < n_x; j++) {
            n = (j - 1);
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        // Internal cells
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

        // Left boundary
        int j = 1;
        for (i = 2; i < n_y; i++) {
            n = (i - 1) * n_x;
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        // Right boundary
        j = n_x;
        for (i = 2; i < n_y; i++) {
            n = i * n_x - 1;
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
            tripletList.emplace_back(n, n + n_x, -A_s(i, j));
        }

        // Bottom boundary
        i = n_y;
        for (j = 2; j < n_x; j++) {
            n = (n_y - 1) * n_x + (j - 1);
            tripletList.emplace_back(n, n, A_p(i, j));
            tripletList.emplace_back(n, n - 1, -A_w(i, j));
            tripletList.emplace_back(n, n + 1, -A_e(i, j));
            tripletList.emplace_back(n, n - n_x, -A_n(i, j));
        }

        // Corners
        // Top-left corner
        n = 0;
        i = 1;
        j = 1;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n + 1, -A_e(i, j));
        tripletList.emplace_back(n, n + n_x, -A_s(i, j));

        // Top-right corner
        i = 1;
        j = n_x;
        n = n_x - 1;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n - 1, -A_w(i, j));
        tripletList.emplace_back(n, n + n_x, -A_s(i, j));

        // Bottom-left corner
        n = (n_y - 1) * n_x;
        j = 1;
        i = n_y;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n + 1, -A_e(i, j));
        tripletList.emplace_back(n, n - n_x, -A_n(i, j));

        // Bottom-right corner
        n = (n_x * n_y) - 1;
        i = n_y;
        j = n_x;
        tripletList.emplace_back(n, n, A_p(i, j));
        tripletList.emplace_back(n, n - 1, -A_w(i, j));
        tripletList.emplace_back(n, n - n_x, -A_n(i, j));

        // Assemble the sparse matrix
        A.setFromTriplets(tripletList.begin(), tripletList.end());
    }
};

void movement_function(Mesh &mesh,Equation &equ_u,Equation &equ_v)

{   double alpha_uv=10e-2;
    int n,i,j;
    int n_x=equ_u.n_x;
    int n_y=equ_u.n_y;
    double D_e,D_w,D_n,D_s,F_e,F_n,F_s,F_w;
    double  re=10000;
    double rho =1;
    D_e=dy/(dx*re);
    D_w=dy/(dx*re);
    D_n=dx/(dy*re);
    D_s=dx/(dy*re);
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
    Eigen::VectorXd &source_x=equ_u.source;
    Eigen::VectorXd &source_y=equ_v.source;
    //内部面
    for(i=2;i<n_y;i++)
    {
        for (j=2;j<n_x;j++)
        {    //nx，每行单元数
            n=(i-1)*n_x + (j-1);
            //计算单元面上流量
            F_e=dy*u_face(i,j);
            F_w=dy*u_face(i,j-1);
            F_n=dx*v_face(i-1,j);
            F_s=dx*v_face(i,j);
            //离散对流项和扩散项 一阶迎风格式
            A_e(i,j)=D_e + std::max(0.0,-F_e);
            A_w(i,j)=D_w + std::max(0.0,F_w);
            A_n(i,j)=D_n + std::max(0.0,-F_n);
            A_s(i,j)=D_s + std::max(0.0,F_s);
            A_p(i,j)=D_e + D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);
            //计算压力梯度，注意矩阵第一行第一列为左上角，亚松弛修正
            source_x[n]=0.5*alpha_uv*(p(i,j-1)-p(i,j+1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);
            source_y[n]=0.5*alpha_uv*(p(i+1,j)-p(i-1,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
        }
    }

    //左壁面
    j=1;
    for (i=2;i<n_y;i++)
    { 
        n=(i-1)*n_x ;

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);  //左壁面速度初始为0
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_s(i,j)=D_s + std::max(0.0,F_s);//面到cell的距离是cell之间的一半
        A_p(i,j)=D_e + 2*D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i - 1,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
    }

    //底部壁面
    i=n_y;
    for(j=2;j<n_x;j++)
    {
        n=(n_y-1)*n_x + (j-1);

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j) ;       

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_p(i,j)=D_e + D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j + 1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);
        source_y[n]=0.5*alpha_uv*(p(i,j) - p(i - 1,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);   
    }
    //右壁面
    j=n_x;
    for(i=2;i<n_y;i++)
    {
        n=i*n_x -1;    
     
        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);   
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_s(i,j)=D_s + std::max(0.0,F_s);
        A_p(i,j)=2*D_e + D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j-1) - p(i,j))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i - 1,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
    }
    //顶盖
    i=1;
    for(j=2;j<n_x;j++)
    {
        n=(j-1);

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_s(i,j)=D_s + std::max(0.0,F_s);
        A_p(i,j)=D_e + D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n) ;

        source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j + 1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx/rho  + (1-alpha_uv)*A_p(i,j)*v_star(i,j);   
    }
    //左上角
    i=1;
    j=1;
    n=0;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1);
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_e(i,j)=D_e + std::max(0.0,-F_e);
    A_s(i,j)=D_s + std::max(0.0,F_s);
    A_p(i,j)=D_e + 2*D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n) ;

    source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));  //P_o - 0.5(P_o+P_e)
    source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);  

    //左下角
    i=1;
    j=n_x;
    n=n_x-1;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1);  
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_w(i,j)=D_w + std::max(0.0,F_w);
    A_s(i,j)=D_s + std::max(0.0,F_s);
    A_p(i,j)=2*D_e + D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));  //0.5(P_w+P_o)-P_o
    source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j); 

    //左下角
    i=n_y;
    j=1;
    n=(n_y-1)*n_x ;


    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1) ; 
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_e(i,j)=D_e + std::max(0.0,-F_e);
    A_n(i,j)=D_n + std::max(0.0,-F_n);
    A_p(i,j)=D_e + 2*D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  //P_o - 0.5(P_o+P_e)
    source_y[n]=0.5*alpha_uv*(p(i,j) - p(i - 1,j))*dx/rho + (1-alpha_uv)*A_p(i,j)*v_star(i,j);  //P_o - 0.5(P_o+P_n)

    //右下角
    i=n_y;
    j=n_x;
    n=n_x*n_y-1;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1) ; 
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_w(i,j)=D_w + std::max(0.0,F_w);
    A_n(i,j)=D_n + std::max(0.0,-F_n);
    A_p(i,j)=2*D_e + D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j - 1)- p(i,j))*dy/rho + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  
    source_y[n]=0.5*alpha_uv*(p(i,j) - p(i-1,j))*dx/rho  + (1-alpha_uv)*A_p(i,j)*v_star(i,j); 


    A_e=alpha_uv*A_e;
    A_w=alpha_uv*A_w;
    A_n=alpha_uv*A_n;
    A_s=alpha_uv*A_s;
    equ_v.A_p=equ_u.A_p;
    equ_v.A_w=equ_u.A_w;
    equ_v.A_e=equ_u.A_e;
    equ_v.A_n=equ_u.A_n;
    equ_v.A_s=equ_u.A_s;


}
//组装方程组


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
    //std::cout << "L2 Norm : " <<  (A*x-b).norm()     << std::endl;
    //std::cout << "Number of iterations : " <<  solver.iterations()    << std::endl;

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
    //u face velocity
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
    
    double alpha_p=0.6;
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


void post_processing(Mesh &mseh,int n_x,int n_y)
{   
    VectorXd x(n_x+2),y(n_y+2);
    x << 0,VectorXd::LinSpaced(n_x,dx/2.0,1-dx/2.0),1;
    y << 0,VectorXd::LinSpaced(n_y,dy/2.0,1-dy/2.0),1;

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


int main ()
{


    std::cout << "x方向上划分个数:";
    std::cin >> n_x0;
    
    std::cout << "y方向上划分个数:";
    std::cin >> n_y0;
    //a，b为网格边长

    std::cout << "正方形网格长度:";
    std::cin >> a;
    
    dx=a/n_x0;
    dy=a/n_y0;
    
    std::cout << "顶盖速度:";
    std::cin >> vx;
   
    velocity=vx;
    //初始化速度压力矩阵
    Mesh mesh(n_y0,n_x0);
    Equation equ_u(n_y0,n_x0);
    Equation equ_v(n_y0,n_x0);
    Equation equ_p(n_y0,n_x0);
   
    
    //设置边界条件   
    mesh.u.block(0,1,1,n_x0)=MatrixXd::Constant(1,n_x0,vx);
    mesh.u_star.block(0,1,1,n_x0)=MatrixXd::Constant(1,n_x0,vx);
    mesh.u_face.block(0,1,1,n_x0-1)=MatrixXd::Constant(1,n_x0-1,vx);
    

    int max_outer_iterations=5000;
    
    for(int n=1;n<=max_outer_iterations;n++)
    {
        movement_function(mesh,equ_u,equ_v);
        equ_u.build_matrix();
        equ_v.build_matrix();
        double epsilon_uv=0.75;
        solve(equ_u.A,equ_u.source,mesh.u,l2_norm_x,epsilon_uv,n_x0,n_y0);
        solve(equ_v.A,equ_v.source,mesh.v,l2_norm_y,epsilon_uv,n_x0,n_y0);
        
        face_velocity(mesh ,equ_u);
        
        pressure_function(mesh,equ_p,equ_u);

        equ_p.build_matrix();
        double epsilon_p=1e-4;
        solve(equ_p.A,equ_p.source,mesh.p_prime,l2_norm_p,epsilon_p,n_x0,n_y0);
        
        correct_pressure(mesh,equ_u);

        correct_velocity(mesh,equ_u);
        
      
        mesh.p=mesh.p_star;
        
        std::cout << " 轮数"<< n << " x速度残差 " << l2_norm_x << " y速度残差 " << l2_norm_y << " 压力残差 " << l2_norm_p << "\n";
        
        if(l2_norm_x < 1e-5 & l2_norm_y < 1e-5 & l2_norm_p < 1e-6)
        { 
        std::cout << "simple算法收敛";
        break;
        }

    }

    post_processing(mesh,n_x0,n_y0);

}