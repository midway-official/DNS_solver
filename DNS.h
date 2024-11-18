// dns.h
#ifndef DNS_H
#define DNS_H

#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <algorithm>
#include <fstream>
#include <vector>

using namespace Eigen;

// 全局变量声明
extern int n_x0, n_y0;
extern double dx, dy, vx;
extern double velocity;
extern double l2_norm_x, l2_norm_y, l2_norm_p;
extern double a, b;

// Mesh 类声明
class Mesh {
public:
    Eigen::MatrixXd u, u_star;
    Eigen::MatrixXd v, v_star;
    Eigen::MatrixXd p, p_star, p_prime;
    Eigen::MatrixXd u_face, v_face;

    Mesh() = default;  // 默认构造函数
    Mesh(int n_y, int n_x);  // 参数化构造函数

    void initializeToZero();
    void displayMatrix(const Eigen::MatrixXd& matrix, const std::string& name) const;
    void displayAll() const;
};

// Equation 类声明
class Equation {
public:
    Eigen::MatrixXd A_p, A_e, A_w, A_n, A_s;
    Eigen::VectorXd source;
    Eigen::SparseMatrix<double> A;
    int n_x, n_y;

    Equation(int n_y_, int n_x_);
    void initializeToZero();
    void build_matrix();
};
void solve(SparseMatrix <double> &A,VectorXd &b,MatrixXd &phi,double &l2_norm,double &epsilon,int n_x,int n_y);

void face_velocity(Mesh &mesh,Equation &equ_u);


void pressure_function(Mesh &mesh,Equation &equ_p,Equation &equ_u);

//修正压力
void correct_pressure(Mesh &mesh,Equation &equ_u);

void correct_velocity(Mesh &mesh,Equation &equ_u);

void post_processing(Mesh &mseh,int n_x,int n_y);

#endif // DNS_H