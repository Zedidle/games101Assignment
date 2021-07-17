//
// Created by LEI XU on 4/11/19.
//

#ifndef RASTERIZER_TRIANGLE_H
#define RASTERIZER_TRIANGLE_H

#include <eigen3/Eigen/Eigen>
#include <vector>
using namespace std;
using namespace Eigen;
class Triangle{
    public:
        Vector3f v[3]; /*the original coordinates of the triangle, v0, v1, v2 in counter clockwise order*/
        /*Per vertex values*/
        Vector3f color[3]; //color at each vertex;
        Vector2f tex_coords[3]; //texture u,v
        Vector3f normal[3]; //normal vector for each vertex

        //Texture *tex;
        Triangle();

        Eigen::Vector3f a() const { return v[0]; }
        Eigen::Vector3f b() const { return v[1]; }
        Eigen::Vector3f c() const { return v[2]; }

        void setVertex(int ind, Vector3f ver); /*set i-th vertex coordinates */
        void setNormal(int ind, Vector3f n); /*set i-th vertex normal vector*/
        void setColor(int ind, float r, float g, float b); /*set i-th vertex color*/
        
        Vector3f getColor() const { return color[0]*255; } // Only one color per triangle.
        Vector3f getColor(float alpha, float beta, float gamma) const { 
            return (color[0]*alpha + color[1]*beta + color[2]*gamma) * 255;
        }
        
        void setTexCoord(int ind, float s, float t); /*set i-th vertex texture coordinate*/
        std::array<Vector4f, 3> toVector4() const;

        vector<Eigen::Vector3f>  getAABB() const{
            float minX, maxX, minY, maxY;
            minX = min(a()[0], min(b()[0], c()[0]));
            maxX = max(a()[0], max(b()[0], c()[0]));
            minY = min(a()[1], min(b()[1], c()[1]));
            maxY = max(a()[1], max(b()[1], c()[1]));
            vector<Eigen::Vector3f> aabb(4);
            aabb[0] = Eigen::Vector3f(minX, minY, 1);
            aabb[1] = Eigen::Vector3f(minX, maxY, 1);
            aabb[2] = Eigen::Vector3f(maxX, maxY, 1);
            aabb[3] = Eigen::Vector3f(maxX, minY, 1);
            return aabb;
        };
};






#endif //RASTERIZER_TRIANGLE_H
