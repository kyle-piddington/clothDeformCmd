#ifndef __MATERIAL_H__
#define __MATERIAL_H__
#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#ifdef __unix__
#include <GL/glut.h>
#endif
#include "GLSL.h"
#define GLM_FORCE_RADIANS

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"


class Material
{
public:
    Material(glm::vec3 ambient, glm::vec3 diffuse, glm::vec3 specular, float specIntensity);
    ~Material();

    /**
     * Bind a material to the current shader
     * @param h_ka Ambient Handle
     * @param h_kd Diffuse Handle
     * @param h_ks Specular Handle
     * @param h_s  Specular Intensity Handle
     */
    void bindMaterial(int h_ka, int h_kd, int h_ks, int h_s);

private:
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float specIntensity;


};

#endif
