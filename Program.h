#pragma  once
#ifndef __Program__
#define __Program__

#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#ifdef __unix__
#include <GL/glut.h>
#endif
#ifdef _WIN32
#define GLFW_INCLUDE_GLCOREARB
#include <GL/glew.h>
#include <cstdlib>
#include <glut.h>
#endif
#include <map>
#include <string>

class Program
{
public:
	Program();
	virtual ~Program();
	
	void setVerbose(bool v) { verbose = v; }
	bool isVerbose() const { return verbose; }
	
	void setShaderNames(const std::string &v, const std::string &f);
	virtual bool init();
	virtual void bind();
	virtual void unbind();

	void addAttribute(const std::string &name);
	void addUniform(const std::string &name);
	GLint getAttribute(const std::string &name) const;
	GLint getUniform(const std::string &name) const;
	
protected:
	std::string vShaderName;
	std::string fShaderName;
	
private:
	GLuint pid;
	std::map<std::string,GLint> attributes;
	std::map<std::string,GLint> uniforms;
	bool verbose;
};

#endif
