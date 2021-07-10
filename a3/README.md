前6项都能得分，最后一项没做

各个函数所实现的功能：
main.cpp Eigen::Matrix4f get_projection_matrix
    正常的透视投影，对正交投影矩阵 M_ortho_zoom 的最后一行乘以-1，达到投影效果 x, y, z都翻转；

main.cpp Eigen::Vector3f phong_fragment_shader
    通过Blinn-phong公式求得环境光，以及每个光源导致的漫反射，高光，所有的光加在一起得到最终该像素的颜色；
    假设光源都是平行光，达到观察点之前不会衰减，因此r_sqr求的观察点离相机距离的平方；
    kd.cwiseProduct(I / r_sqr) 的cwiseProduct用的就很巧妙；
    最后调整光源l1, l2达到目标效果；

main.cpp Eigen::Vector3f texture_fragment_shader
    同理于phong_fragment_shader，另外只需要将payload中的texture以及纹理坐标求得对应纹理位置的颜色值并存储在return_color，然后替换漫反射系数kd；


main.cpp displacement_fragment_shader
    通过公式求出TBN；
    texture->getGrayColor：由于没有高度图而套用灰度公式；
    最终求得观察点的位置偏移以及改变后的法线来作Blinn-phong计算；

main.cpp bump_fragment_shader
    同理于displacement_fragment_shader，只是没有观察点并没有实际的位移；


