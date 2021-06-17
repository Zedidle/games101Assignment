#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

using namespace std;
using namespace Eigen;

#define PI acos(-1)

// myutil
namespace mu{

    float radian(float value){
        return value/180 * PI;
    }

    float angleVector3f(Vector3f a, Vector3f b){
        float cosValue = a.dot(b) / (a.norm() * b.norm());
        return acos(cosValue) / PI * 180;
    }

    Vector3f projection(Vector3f a, Vector3f b){
        float angle = angleVector3f(a, b);
        float bLength = cos(angle/180 * PI) * a.norm();
        return bLength * b.normalized();
    }

};