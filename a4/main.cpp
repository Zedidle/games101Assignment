#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    // if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 4) 
    // {
    //     std::cout << "Left button of the mouse is clicked - position (" << x << ", "
    //     << y << ")" << '\n';
    //     control_points.emplace_back(x, y);
    // }
    if (event == cv::EVENT_LBUTTONDOWN){
        control_points.emplace_back(x, y);
    }
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &points, float t) 
{
    // TODO: Implement de Casteljau's algorithm
    if(points.size() < 2) return cv::Point2f(0,0); 
    if(t > 1) return points[1];
    std::vector<cv::Point2f> P = points;

    while (P.size()>2){
        for(int i=0;i<P.size()-1;i++){
            P[i] = P[i] + t * (P[i+1] - P[i]);
        }
        P.pop_back();
    }
    return P[0] + t * (P[1] - P[0]);
}

void bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        cv::Point2f point = recursive_bezier(points, t);
        window.at<cv::Vec3b>(point.y, point.x)[4] = 255;
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);
    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 1, {255, 255, 255}, 3);
        }

        bezier(control_points, window);
        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}
