#include "Light.h"


Light::Light(glm::vec3 position, float intensity):
    position(position),
    intensity(intensity){}

Light::~Light(){}


void Light::setPosition(glm::vec3 newPos){
    this->position = newPos;
}
void Light::movePosition(glm::vec3 delt){
    this->position += delt;
}

