#include <iostream>
#include <opencv2/opencv.hpp>

#include "global.hpp"
#include "rasterizer.hpp"
#include "Triangle.hpp"
#include "Shader.hpp"
#include "Texture.hpp"
#include "OBJ_Loader.h"
#include "../mylib.h"
#include "math.h"

Eigen::Vector3f eye_pos = {0,0,10};
using namespace std;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float angle)
{
    Eigen::Matrix4f rotation;
    angle = angle * MY_PI / 180.f;
    rotation << cos(angle), 0, sin(angle), 0,
                0, 1, 0, 0,
                -sin(angle), 0, cos(angle), 0,
                0, 0, 0, 1;

    Eigen::Matrix4f scale;
    scale << 2.5, 0, 0, 0,
              0, 2.5, 0, 0,
              0, 0, 2.5, 0,
              0, 0, 0, 1;

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    return translate * rotation * scale;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // TODO: Use the same projection matrix from the previous assignments
    float Tn = tan(mu::radian(eye_fov/2)) * zNear;
    float Tf = tan(mu::radian(eye_fov/2)) * zFar;
    float Rn = aspect_ratio * Tn;
    float Rf = aspect_ratio * Tf;

    Matrix4f M_persp2ortho;
    M_persp2ortho << zNear,0,0,0,  0,zNear,0,0,  0,0,zNear+zFar,-zNear*zFar,  0,0,1,0;

    Matrix4f M_ortho, M_ortho_zoom, M_ortho_translation; // Zoom, Translation
    M_ortho_zoom << 1/Rn,0,0,0, 0,1/Tn,0,0, 0,0,1/zNear,0, 0,0,0,-1/3;
    M_ortho_translation << 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1;
    M_ortho = M_ortho_zoom * M_ortho_translation;

    return M_ortho * M_persp2ortho;
}

Eigen::Vector3f vertex_shader(const vertex_shader_payload& payload)
{
    return payload.position;
}

Eigen::Vector3f normal_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f return_color = (payload.normal.head<3>().normalized() + Eigen::Vector3f(1.0f, 1.0f, 1.0f)) / 2.f;
    Eigen::Vector3f result;
    result << return_color.x() * 255, return_color.y() * 255, return_color.z() * 255;
    return result;
}

static Eigen::Vector3f reflect(const Eigen::Vector3f& vec, const Eigen::Vector3f& axis)
{
    auto costheta = vec.dot(axis);
    return (2 * costheta * axis - vec).normalized();
}

struct light
{
    Eigen::Vector3f position;
    Eigen::Vector3f intensity;
};

Eigen::Vector3f texture_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f return_color = {0, 0, 0};
    if (payload.texture)
    {
        // TODO: Get the texture value at the texture coordinates of the current fragment
        return_color = payload.texture->getColor(payload.tex_coords[0], payload.tex_coords[1]);
    }

    Eigen::Vector3f texture_color = return_color;

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = texture_color / 255.f;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {300, 300, 300}};
    auto l2 = light{{-20, 20, 0}, {200, 200, 200}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};

    float p = 150;

    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};

    for (auto& light : lights)
    {
        Eigen::Vector3f l = light.position - point;  // PointToLight 
        Eigen::Vector3f v = eye_pos - point;  // PointToEye
        float r_sqr = v.dot(v);
        v.normalize();
        l.normalize();

        Eigen::Vector3f I = light.intensity;   // 光强
        Eigen::Vector3f h = (l+v).normalized();  // 半程向量

        Eigen::Vector3f Ld = kd.cwiseProduct(I / r_sqr) * std::max(0.0f, normal.dot(l));
        Eigen::Vector3f Ls = ks.cwiseProduct(I / r_sqr) * std::pow(std::max(0.0f, normal.dot(h)), p);
        result_color += Ld + Ls;
    }
    Eigen::Vector3f La = ka.cwiseProduct(amb_light_intensity);

    return result_color * 255.f;
}


Eigen::Vector3f kd_Falloff = Eigen::Vector3f(0,0,0);
Eigen::Vector3f phong_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f kd = payload.color;                             // 颜色 - 漫反射系数
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);   // 高光系数
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);      // 环境光系数

    auto l1 = light{{20, 20, 20}, {100, 100, 100}};
    auto l2 = light{{-20, 20, 0}, {200, 200, 200}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};

    float p = 64; // 高光指数

    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};
    for (auto& light : lights)
    {
        Eigen::Vector3f l = light.position - point;  // PointToLight 
        Eigen::Vector3f v = eye_pos - point;  // PointToEye
        float r_sqr = v.dot(v);
        v.normalize();
        l.normalize();

        Eigen::Vector3f I = light.intensity;   // 光强
        Eigen::Vector3f h = (l+v).normalized();  // 半程向量

        Eigen::Vector3f Ld = kd.cwiseProduct(I / r_sqr) * std::max(0.0f, normal.dot(l));
        Eigen::Vector3f Ls = ks.cwiseProduct(I / r_sqr) * std::pow(std::max(0.0f, normal.dot(h)), p);
        result_color += Ld + Ls;
    }
    Eigen::Vector3f La = ka.cwiseProduct(amb_light_intensity);

    return (result_color + La) * 255.f;
}


Eigen::Vector3f displacement_fragment_shader(const fragment_shader_payload& payload)
{
    
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {200, 200, 200}};
    auto l2 = light{{-20, 20, 0}, {50, 50, 50}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color; 
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;

    Eigen::Vector3f n = normal;
    float x = n[0],  y = n[1], z = n[2];
    float sqxz = sqrt(x*x+z*z);
    Eigen::Vector3f t =  Eigen::Vector3f(x*y/sqxz, sqxz, z*y/sqxz);
    t = (t - n.dot(t) * n).normalized();
    Eigen::Vector3f b = n.cross(t).normalized();
    Eigen::Matrix3f TBN;
    TBN << t, b, n; // TBN的定义？

    Texture* texture = payload.texture;
    float w = texture->width;
    float h = texture->height;
    float u = payload.tex_coords[0];
    float v = payload.tex_coords[1];

    float dU = kh * kn * (texture->getGrayColor(u+1/w,v) - texture->getGrayColor(u,v));
    float dV = kh * kn * (texture->getGrayColor(u,v+1/h) - texture->getGrayColor(u,v));
    Eigen::Vector3f ln = {-dU, -dV, 1};
    point = point + kn * n * texture->getGrayColor(u,v);
    normal = (TBN * ln).normalized();

    Eigen::Vector3f result_color = {0, 0, 0};

    for (auto& light : lights)
    {
        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
        Eigen::Vector3f l = light.position - point;  // PointToLight 
        Eigen::Vector3f v = eye_pos - point;  // PointToEye
        float square_r = v.dot(v);
        l.normalize();
        v.normalize();

        Eigen::Vector3f I = light.intensity;   // 光强
        Eigen::Vector3f h = (l + v).normalized();  // 半程向量

        Eigen::Vector3f Ld = kd.cwiseProduct(I / square_r) * std::max(0.0f, normal.dot(l));
        Eigen::Vector3f Ls = ks.cwiseProduct(I / square_r) * std::pow(std::max(0.0f, normal.dot(h)), p);
        
        result_color += Ld + Ls;
    }

    Eigen::Vector3f La = ka.cwiseProduct(amb_light_intensity);
    return (result_color + La) * 255.f;
}


Eigen::Vector3f bump_fragment_shader(const fragment_shader_payload& payload)
{
    
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color; 
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;  // kh，kn又是什么？

    Eigen::Vector3f n = normal;
    float x = n[0],  y = n[1], z = n[2];
    float sqxz = sqrt(x*x+z*z);
    Eigen::Vector3f t =  Eigen::Vector3f(x*y/sqxz, sqxz, z*y/sqxz);
    t = (t - n.dot(t) * n).normalized();
    Eigen::Vector3f b = n.cross(t).normalized();
    Eigen::Matrix3f TBN;
    TBN << t, b, n; // TBN的定义？

    Texture* texture = payload.texture;
    float w = texture->width;
    float h = texture->height;
    float u = payload.tex_coords[0]; 
    float v = payload.tex_coords[1];

    float dU = kh * kn * (texture->getGrayColor(u+1/w, v) - texture->getGrayColor(u,v));
    float dV = kh * kn * (texture->getGrayColor(u, v+1/h) - texture->getGrayColor(u,v));
    Eigen::Vector3f ln = {-dU, -dV, 1};
    normal = (TBN * ln).normalized();

    Eigen::Vector3f result_color = {0, 0, 0};
    result_color = normal;
    return result_color * 255.f;
}

int main(int argc, const char** argv)
{
    std::vector<Triangle*> TriangleList;

    float angle = 140.0;
    bool command_line = false;
    bool shader_all = false;

    std::string filename = "output.png";
    string shader_type;
    objl::Loader Loader;
    std::string obj_path = "../models/spot/";

    // Load .obj File
    vector<string> objs = {
        "../models/spot/spot_triangulated_good.obj",
        "../models/cube/cube.obj",
        "../models/Crate/Crate1.obj",
        "../models/bunny/bunny.obj"
    };
    int obj_index = 3;
    bool loadout = Loader.LoadFile(objs[obj_index]);
    for(auto mesh:Loader.LoadedMeshes)
    {
        for(int i=0;i<mesh.Vertices.size();i+=3)
        {
            Triangle* t = new Triangle();
            for(int j=0;j<3;j++)
            {
                t->setVertex(j,Vector4f(mesh.Vertices[i+j].Position.X,mesh.Vertices[i+j].Position.Y,mesh.Vertices[i+j].Position.Z,1.0));
                t->setNormal(j,Vector3f(mesh.Vertices[i+j].Normal.X,mesh.Vertices[i+j].Normal.Y,mesh.Vertices[i+j].Normal.Z));
                t->setTexCoord(j,Vector2f(mesh.Vertices[i+j].TextureCoordinate.X, mesh.Vertices[i+j].TextureCoordinate.Y));
            }
            TriangleList.push_back(t);
        }
    }

    rst::rasterizer r(700, 700);

    auto texture_path = "hmap.jpg";
    r.set_texture(Texture(obj_path + texture_path));

    std::function<Eigen::Vector3f(fragment_shader_payload)> active_shader = phong_fragment_shader;
    map<string, function<Eigen::Vector3f(fragment_shader_payload)>> shader_Type2Func = {
        {"phong", phong_fragment_shader},
        {"texture", texture_fragment_shader},
        {"normal", normal_fragment_shader},
        {"bump", bump_fragment_shader},
        {"displacement", displacement_fragment_shader}
    };

    if (argc >= 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
        shader_type = std::string(argv[2]);

        if (argc == 3 && shader_type == "all")
        {
            shader_all = true;
        }
        else if(argc==3){
            for(auto& p: shader_Type2Func){
                if(p.first == "texture"){
                    texture_path = "spot_texture.png";
                    r.set_texture(Texture(obj_path + texture_path));
                }
                if(p.first == shader_type){
                    active_shader = p.second;
                    break;
                }
            }
        }
    }

    int key = 0;
    int frame_count = 0;


    r.set_vertex_shader(vertex_shader);
    r.set_fragment_shader(active_shader);

    if (command_line)
    {
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        if(shader_all){
            for(auto p : shader_Type2Func){
                if(p.first == "texture"){
                    texture_path = "spot_texture.png";
                    r.set_texture(Texture(obj_path + texture_path));
                }

                r.set_fragment_shader(p.second);
                r.draw(TriangleList);
                cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
                image.convertTo(image, CV_8UC3, 1.0f);
                cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
                cv::imwrite(filename+"/"+p.first+".jpg", image);
                r.clear(rst::Buffers::Color | rst::Buffers::Depth);
            }
        }else{
            

            r.draw(TriangleList);
            cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
            image.convertTo(image, CV_8UC3, 1.0f);
            cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
            cv::imwrite(filename, image); 
        }

        return 0;
    }

    while(key != 27)
    {
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);
        cv::imwrite(filename, image);
        key = cv::waitKey(10);
        
        // cout << "frame_count:" << frame_count++ << endl; 
        // angle += 5;
        if (key == 'a' )
        {
            angle -= 0.1;
        }
        else if (key == 'd')
        {
            angle += 0.1;
        }

    }
    return 0;
}
