//
// Created by goksu on 4/6/19.
//

#pragma once

#include "shapes/Triangle.hpp"
#include "shapes/Rectangle.hpp"
#include <algorithm>
#include <eigen3/Eigen/Eigen>
using namespace Eigen;
using namespace std;

namespace rst {
    enum class Buffers
    {
        Color = 1,
        Depth = 2
    };

    inline Buffers operator|(Buffers a, Buffers b)
    {
        return Buffers((int)a | (int)b);
    }

    inline Buffers operator&(Buffers a, Buffers b)
    {
        return Buffers((int)a & (int)b);
    }

    enum class Primitive
    {
        Line,
        Triangle,
        Rectangle,
    };

    /*
    * For the curious : The draw function takes two buffer id's as its arguments.
    * These two structs make sure that if you mix up with their orders, the
    * compiler won't compile it. Aka : Type safety
    */
    struct pos_buf_id
    {
        int pos_id = 0;
    };

    struct ind_buf_id
    {
        int ind_id = 0;
    };

    class rasterizer
    {
        public:
            rasterizer(int w, int h);

            pos_buf_id load_positions(const vector<Eigen::Vector3f>& positions);
            ind_buf_id load_indices(const vector<Eigen::Vector3i>& indices);

            void set_model(const Eigen::Matrix4f& m);
            void set_view(const Eigen::Matrix4f& v);
            void set_projection(const Eigen::Matrix4f& p);

            void clear(Buffers buff);

            void draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, Primitive type);
            void drawTriangle(vector<Eigen::Vector3f>& buf, vector<Eigen::Vector3i>& ind);
            void drawRectangle(vector<Eigen::Vector3f>& buf, vector<Eigen::Vector3i>& ind);

            vector<Eigen::Vector3f>& frame_buffer() { return frame_buf; }
            
            void draw_line(Eigen::Vector3f begin, Eigen::Vector3f end);

        private:
            void set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color);

        private:
            Eigen::Matrix4f model;
            Eigen::Matrix4f view;
            Eigen::Matrix4f projection;

            map<int, vector<Eigen::Vector3f>> pos_buf;
            map<int, vector<Eigen::Vector3i>> ind_buf;

            vector<Eigen::Vector3f> frame_buf;
            vector<float> depth_buf;
            int get_index(int x, int y);

            int width, height;

            int next_id = 0;
            int get_next_id() { return next_id++; }
    };
} // namespace rst