#include "DNS.h"
#include <filesystem>
#include <chrono>
namespace fs = std::filesystem;


void show_progress_bar(int current_step, int total_steps, double elapsed_time) {
    // 计算进度百分比
    double progress = static_cast<double>(current_step) / total_steps;
    
    // 设置进度条的宽度
    int bar_width = 50;
    
    // 计算进度条中"="的数量
    int pos = static_cast<int>(bar_width * progress);
    
    // 计算预计剩余时间
    double remaining_time = (elapsed_time / current_step) * (total_steps - current_step);
    
    // 打印进度条和相关信息
    std::cout << "[";
    
    // 绘制进度条
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) {
            std::cout << "=";  // 已完成的部分
        } else if (i == pos) {
            std::cout << ">";  // 当前进度的位置
        } else {
            std::cout << " ";  // 未完成的部分
        }
    }
    
    // 显示进度条，已用时间和预计剩余时间
    std::cout << "] " 
              << std::fixed << std::setprecision(2) << progress * 100 << "% "  // 显示进度百分比
              << "已用时间: " << std::fixed << std::setprecision(2) << elapsed_time << "秒 "  // 显示已用时间
              << "预计剩余时间: " << std::fixed << std::setprecision(2) << remaining_time << "秒\r";  // 显示预计剩余时间
    
    // 刷新输出，确保实时更新
    std::cout.flush();
}
void movement_function(Mesh &mesh,Equation &equ_u,Equation &equ_v)

{   double alpha_uv=10e-2;
    int n,i,j;
    int n_x=equ_u.n_x;
    int n_y=equ_u.n_y;
    double D_e,D_w,D_n,D_s,F_e,F_n,F_s,F_w;
    double  vp=10e-6;
    double rho =10e3;
    D_e=(dy*vp)/(dx);
    D_w=(dy*vp)/(dx);
    D_n=(dx*vp)/(dy);
    D_s=(dx*vp)/(dy);
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
    //时间步长
    double dt;
    std::cout << "时间步长:";
    std::cin >> dt;
    //总时间步数
    int timesteps;
    std::cout << "时间步长数:";
    std::cin >> timesteps;
    //设定全局变量
    velocity=vx;
    //创建网格
    Mesh mesh(n_y0,n_x0);
    //建立u v p的方程
    Equation equ_u(n_y0,n_x0);
    Equation equ_v(n_y0,n_x0);
    Equation equ_p(n_y0,n_x0);
    
    
    //设置边界条件   
    mesh.u.block(0,1,1,n_x0)=MatrixXd::Constant(1,n_x0,vx);
    mesh.u_star.block(0,1,1,n_x0)=MatrixXd::Constant(1,n_x0,vx);
    mesh.u_face.block(0,1,1,n_x0-1)=MatrixXd::Constant(1,n_x0-1,vx);
    
        // 创建顶层的 result 文件夹
    std::string result_folder = "result";
    if (!fs::exists(result_folder)) {
            fs::create_directory(result_folder);
            
        }

    // 切换到 result 文件夹
    fs::current_path(result_folder);
    auto start_time = std::chrono::steady_clock::now();  // 开始计时
    // 循环执行
    for (int i = 0; i <= timesteps; ++i) { 
        std::string folder_name = std::to_string(i);
        if (!fs::exists(folder_name)) {
            fs::create_directory(folder_name);
            std::cout << "结果保存于: " << folder_name << std::endl;
        }
        std::cout<<"时间步长 "<< i <<std::endl;
        // 切换到当前编号文件夹
       fs::current_path(folder_name);
       //记录上一个时间步长的u v
       mesh.u0 = mesh.u;
       mesh.v0 = mesh.v;
       //每步最大迭代次数
       int max_outer_iterations=300;
       //simple算法迭代
       for(int n=1;n<=max_outer_iterations;n++)
          {
          //离散动量方程 
          movement_function(mesh,equ_u,equ_v);
          //组合线性方程组
          equ_u.build_matrix();
          equ_v.build_matrix();
          //求解线性方程组
          double epsilon_uv=0.75;
          solve(equ_u.A,equ_u.source,mesh.u,l2_norm_x,epsilon_uv,n_x0,n_y0);
          solve(equ_v.A,equ_v.source,mesh.v,l2_norm_y,epsilon_uv,n_x0,n_y0);
          //速度插值到面
          face_velocity(mesh ,equ_u);
          //离散压力修正方程
          pressure_function(mesh,equ_p,equ_u);
          //组合线性方程组
          equ_p.build_matrix();
          //求解压力修正方程
          double epsilon_p=1e-4;
          solve(equ_p.A,equ_p.source,mesh.p_prime,l2_norm_p,epsilon_p,n_x0,n_y0);
          //压力修正
          correct_pressure(mesh,equ_u);
          //速度修正
          correct_velocity(mesh,equ_u);
        
          //更新压力
          mesh.p=mesh.p_star;
          //收敛性判断
           std::cout << std::scientific 
          << " 轮数 " << n 
          << " x速度残差 " << std::setprecision(6) << l2_norm_x 
          << " y速度残差 " << std::setprecision(6) << l2_norm_y 
          << " 压力残差 " << std::setprecision(6) << l2_norm_p 
          << "\n" << std::endl;
           if(l2_norm_x < 1e-5 & l2_norm_y < 1e-5 & l2_norm_p < 1e-6)
           { 
           //std::cout << "simple算法收敛"<<std::endl;
            break;
           }
        
        }
    //显式时间推进
    mesh.u = mesh.u0 + dt*mesh.u;
    mesh.v = mesh.v0 + dt*mesh.v;
    // 显示进度条
    auto elapsed_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    show_progress_bar(i, timesteps, elapsed_time);
    //保存信息
    post_processing(mesh,n_x0,n_y0,a);
     // 返回到 result 文件夹
    fs::current_path("..");
    }
        // 返回到程序主目录
    fs::current_path("..");
    // 最后显示实际计算总耗时
    auto total_elapsed_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
    std::cout << "\n计算完成 总耗时: " << total_elapsed_time << "秒" << std::endl;

}