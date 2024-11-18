#include "DNS.h"



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