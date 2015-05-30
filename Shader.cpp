#include "Shader.h"

Shader::Shader(std::string name, int pid, GLenum shaderType){
    GLint rc;
    sHandler = glCreateShader(shaderType);
    const char *shader = GLSL::textFileRead(name.c_str());
    glShaderSource(sHandler,1,&shader,NULL);
    glCompileShader(sHandler);
    GLSL::printError();
    glGetShaderiv(sHandler, GL_COMPILE_STATUS, &rc);
    GLSL::printShaderInfoLog(sHandler);
    if(!rc) {
        printf("Error compiling shader %s\n", name.c_str());
        printf("GL Version: %s\n", glGetString(GL_VERSION));
    }
    this->pid = pid;
}
Shader::~Shader(){}

GLint Shader::bindAttribute(const char * name){
    return GLSL::getAttribLocation(pid, name);
}

GLint Shader::bindUniform(const char* name){
    return GLSL::getUniformLocation(pid, name);
}

void Shader::attach(){
    glAttachShader(pid, this->sHandler);
}
void Shader::detatch(){
    glDetachShader(pid, this->sHandler);
}
