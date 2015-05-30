#include "Material.h"

Material::Material(glm::vec3 ambient, glm::vec3 diffuse, glm::vec3 specular, float specIntensity):
    ambient(ambient),
    diffuse(diffuse),
    specular(specular),
    specIntensity(specIntensity)
    {}

Material::~Material(){}

void Material::bindMaterial(int h_ka, int h_kd, int h_ks, int h_s){
    glUniform3fv(h_ka,  1, glm::value_ptr(ambient));
    glUniform3fv(h_kd,  1, glm::value_ptr(diffuse));
    glUniform3fv(h_ks, 1, glm::value_ptr(specular));
    glUniform1f(h_s, specIntensity);
}
