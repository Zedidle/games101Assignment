#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

#include "../mylib.h"

#define PI acos(-1)

using namespace std;
using namespace Eigen;



void ExampleCpp();
void ExampleVector();
void ExampleMatrix();



int main(){

    // ExampleCpp();
    // ExampleVector();
    ExampleMatrix();

    return 0;
}

void ExampleCpp(){
    cout << "Example of cpp \n";
    // cout << sqrt(2) << endl;
    cout << PI << endl;
    cout << cos(mu::radian(180)) << endl;
    cout << cos(0) << endl;
    cout << sin(mu::radian(30)) << endl;
    cout << endl;
}

void ExampleVector(){
    cout << "Example of vector \n";
    // vector definition
    // Vector3f v(sqrt(3),0.0f,1.0f);
    // Vector3f v(1.0f, 0, 1.0f);
    // Vector3f v(1.0f,0.0f,1.0f);
    Vector3f v(0.0f,1.0f,-1.0f);
    // Vector3f w(sqrt(3),1.0f,0.0f);
    // Vector3f w(1.0f, sqrt(2) / sqrt(3), 1.0f);
    // Vector3f w(0.0f,1.0f,0.0f);
    Vector3f w(1.0f,-1.0f,0.0f);  
    // cout << v.cross(w) << endl;

    // vector output
    // cout << v << endl;
    // v.normalize();
    cout << w.norm() << endl; // 求长度
    // cout << v.squaredNorm() << endl;
    // cout << v.normalized() << endl;
    // cout << v.transpose() << endl;

    // vector add
    // cout << "Example of add \n";
    // cout << v + w << endl;

    // vector scalar multiply
    // cout << "Example of scalar multiply \n";
    // cout << v * 3.0f << endl;
    // cout << 2.0f * v << endl;
    // cout << v.dot(w) << endl;



}

void ExampleMatrix(){
    cout << "Example of matrix \n";

    Vector3f w(1.0f,1.0f,1.0f);
    // define a triangle
    // Matrix3f m;
    // m << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    // cout << m * w << endl;

    // matrix definition
    Matrix3f i,j;
    i << 1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0;
    j << 2.0,3.0,1.0, 4.0,6.0,5.0, 9.0,7.0,8.0;

    Matrix3f normal2;
    normal2 << 2,0,0, 0,2,0, 0,0,2;
    // matrix output
    cout << "Example of output \n";
    // cout << i << endl;
    // matrix add i + j
    // cout << i + j << endl << endl;


    // matrix scalar multiply i * 2.0
    // cout << i * 2 << endl << endl;
    // cout << i * normal2 << endl << endl;

    // matrix multiply i * j
    // cout << i * j << endl << endl;


    float a = mu::radian(90);
    Matrix4f origin;
    Matrix4f RotationZ;
    origin << 1,1,1,0, 2,2,2,0, 3,3,3,0, 0,0,0,1;
    RotationZ << cos(a),-sin(a),0,0,  sin(a),cos(a),0,0, 0,0,1,0, 0,0,0,1;

    cout << origin << endl;
    cout << RotationZ << endl;
    cout << RotationZ * origin << endl;


    
}