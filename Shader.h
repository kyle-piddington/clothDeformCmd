#ifndef __SHADER_H__
#define __SHADER_H__
#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#ifdef __unix__
#include <GL/glut.h>
#endif
#include "GLSL.h"
#include "glm/glm.hpp"
#include <iostream>
class Shader
{
public:
    Shader(std::string name, int pid, GLenum shaderType);
    ~Shader();

    /**
     * Bind an attribute to the sahder
     * @param  name  Name of the attribute in the GLSL program
     * @return       Handler index;
     */
    /**
     * Attatch an attribute to the shader.
     * @param  name [Name of the attribute]
     * @return      [Shader handler]
     */
    int bindAttribute(const char* name);

    /**
     * Bind a uniform to the shader
     * @param  name Name of the uniform
     * @return      Uniform handle
     */
    int bindUniform(const char* name);

    /**
     * Attatch the shader to the program.
     */
    void attach();

    /**
     * Detatch the shader from the program.
     */
    void detatch();
private:
    int pid;
    GLuint sHandler;


};

#endif
