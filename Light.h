#ifndef __LIGHT_H__
#define __LIGHT_H__
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
class Light

{
public:

    Light(glm::vec3 position, float intensity);
    ~Light();

    /**
     * Set the light position to a direct position
     * @param newPos New position
     */
    void setPosition(glm::vec3 newPos);

    /**
     * Move the position by an offset
     * @param delta A 3 vector offset
     */
    void movePosition(glm::vec3 delta);

    /**
     * Current position of the light in world space.
     */
    glm::vec3 position;

    /**
     * Light intensity
     */
    float intensity;
};
#endif
