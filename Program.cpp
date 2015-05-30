#include <iostream>
#include <cassert>
#include <stdio.h>
#include "Program.h"
#include "GLSL.h"

using namespace std;

Program::Program() :
	pid(0),
	vShaderName(""),
	fShaderName(""),
	verbose(true)
{
	
}

Program::~Program()
{
	
}

void Program::setShaderNames(const string &v, const string &f)
{
	vShaderName = v;
	fShaderName = f;
}

bool Program::init()
{
	GLint rc;
	
	// Create shader handles
	GLuint VS = glCreateShader(GL_VERTEX_SHADER);
	GLuint FS = glCreateShader(GL_FRAGMENT_SHADER);
	
	// Read shader sources
	const char *vshader = GLSL::textFileRead(vShaderName.c_str());
	const char *fshader = GLSL::textFileRead(fShaderName.c_str());
	glShaderSource(VS, 1, &vshader, NULL);
	glShaderSource(FS, 1, &fshader, NULL);
	
	// Compile vertex shader
	glCompileShader(VS);
	glGetShaderiv(VS, GL_COMPILE_STATUS, &rc);
	if(!rc) {
		if(isVerbose()) {
			GLSL::printShaderInfoLog(VS);
			printf("Error compiling vertex shader %s\n", vShaderName.c_str());
		}
		return false;
	}
	
	// Compile fragment shader
	glCompileShader(FS);
	glGetShaderiv(FS, GL_COMPILE_STATUS, &rc);
	if(!rc) {
		if(isVerbose()) {
			GLSL::printShaderInfoLog(FS);
			printf("Error compiling fragment shader %s\n", fShaderName.c_str());
		}
		return false;
	}
	
	// Create the program and link
	pid = glCreateProgram();
	glAttachShader(pid, VS);
	glAttachShader(pid, FS);
	glLinkProgram(pid);
	glGetProgramiv(pid, GL_LINK_STATUS, &rc);
	if(!rc) {
		if(isVerbose()) {
			GLSL::printProgramInfoLog(pid);
			printf("Error linking shaders %s and %s\n", vShaderName.c_str(), fShaderName.c_str());
		}
		return false;
	}
	
	GLSL::printError();
	assert(glGetError() == GL_NO_ERROR);
	return true;
}

void Program::bind()
{
	glUseProgram(pid);
}

void Program::unbind()
{
	glUseProgram(0);
}

void Program::addAttribute(const string &name)
{
	attributes[name] = GLSL::getAttribLocation(pid, name.c_str());
}

void Program::addUniform(const string &name)
{
	uniforms[name] = GLSL::getUniformLocation(pid, name.c_str());
}

GLint Program::getAttribute(const string &name) const
{
	map<string,GLint>::const_iterator attribute = attributes.find(name.c_str());
	if(attribute == attributes.end()) {
		if(isVerbose()) {
			cout << name << " is not an attribute variable" << endl;
		}
		return 0;
	}
	return attribute->second;
}

GLint Program::getUniform(const string &name) const
{
	map<string,GLint>::const_iterator uniform = uniforms.find(name.c_str());
	if(uniform == uniforms.end()) {
		if(isVerbose()) {
			cout << name << " is not a uniform variable" << endl;
		}
		return 0;
	}
	return uniform->second;
}
