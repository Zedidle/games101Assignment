#include "shapes/Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "../mylib.h"

using namespace Eigen;
constexpr double MY_PI = 3.1415926;
cv::Point MOUSEMOVE_DIS = cv::Point(0,0);
Matrix4f M_VIEW_ROTATION = Matrix4f::Identity();

Matrix4f get_view_matrix(Vector3f eye_pos)
{
    float a = mu::radian(-MOUSEMOVE_DIS.y);
    float b = mu::radian(-MOUSEMOVE_DIS.x);
    Matrix4f RotationX, RotationY;
    RotationX << 1,0,0,0,  0,cos(a),-sin(a),0,  0,sin(a),cos(a),0,  0,0,0,1;
    RotationY << cos(b),0,sin(b),0,  0,1,0,0,  -sin(b),0,cos(b),0,  0,0,0,1;
    M_VIEW_ROTATION = RotationX * RotationY * M_VIEW_ROTATION;
    Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],  
                 0,1,0,-eye_pos[1],  
                 0,0,1,-eye_pos[2],  
                 0,0,0,1;
    return M_VIEW_ROTATION * translate;
}

Matrix4f get_model_matrix(float rotation_angle)
{

    Matrix4f model = Matrix4f::Identity();
    Matrix4f RotationX, RotationY, RotationZ;
    float a = mu::radian(rotation_angle);
    RotationX << 1,0,0,0,  0,cos(a),-sin(a),0,  0,sin(a),cos(a),0,  0,0,0,1;
    RotationY << cos(a),0,sin(a),0,  0,1,0,0,  -sin(a),0,cos(a),0,  0,0,0,1;
    RotationZ << cos(a),-sin(a),0,0,  sin(a),cos(a),0,0,  0,0,1,0,  0,0,0,1;

    model = RotationZ;
    return model;
}

Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    float Tn = tan(mu::radian(eye_fov/2)) * zNear;
    float Tf = tan(mu::radian(eye_fov/2)) * zFar;
    float Rn = aspect_ratio * Tn;
    float Rf = aspect_ratio * Tf;

    Matrix4f M_ortho, M_persp2ortho;
    Matrix4f M_ortho_zoom, M_ortho_translation; // Zoom, Translation
    // Why this way is right?
    M_ortho_zoom << 1/Rn,0,0,0, 0,1/Tn,0,0, 0,0,1/zNear,0, 0,0,0,1;
    M_ortho_translation << 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1;
    M_ortho = M_ortho_zoom * M_ortho_translation;

    M_persp2ortho << zNear,0,0,0,  0,zNear,0,0,  0,0,zNear+zFar,-zNear*zFar,  0,0,1,0;
    return M_ortho * M_persp2ortho;
}

cv::Point previousPoint = cv::Point(0,0);
void On_mouse(int event, int x, int y, int flags, void*)
{
    return;
	if (event == cv::EVENT_MOUSEMOVE)
	{
        if(previousPoint != cv::Point(0,0)){
            if(abs(x-previousPoint.x) < 3){
                MOUSEMOVE_DIS.x = x-previousPoint.x;
                cout << "MOUSEMOVE_DIS.x: "<< MOUSEMOVE_DIS.x << endl;
            }
            if(abs(y-previousPoint.y) < 3){
                MOUSEMOVE_DIS.y = y-previousPoint.y;
                cout << "MOUSEMOVE_DIS.y: "<< MOUSEMOVE_DIS.y << endl;
            }
        }
		previousPoint = cv::Point(x, y);
	}
}


int main(int argc, const char** argv)
{

    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);
    Vector3f eye_pos = {0, 0, 5};

    std::vector<Vector3f> pos{{2, 0, 0}, {0, 3, 2}, {-2, 0, 2}};
    std::vector<Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 16/9, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }


    while (key != 27) { // key == 27 => ESC
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        
        // maybe if it's buffers over the window'size, it will break;
        // eye_pos = {5, 5,  5 + sin(mu::radian(6*frame_count)) };

        // float eye_forv = 120 - frame_count % 60;
        float eye_forv = 120;
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(eye_forv, 16/9, 1, 100));
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
	    // cv::setMouseCallback("image", On_mouse, 0);
        key = cv::waitKey(1);

        switch(key){
            case 'a': eye_pos.x() += 0.1; break;
            case 'd': eye_pos.x() -= 0.1; break;
            case 'w': eye_pos.z() -= 0.1; break;
            case 's': eye_pos.z() += 0.1; break;
            case 'q': eye_pos.y() += 0.1; break;
            case 'e': eye_pos.y() -= 0.1; break;
        }
    }

    return 0;
}
