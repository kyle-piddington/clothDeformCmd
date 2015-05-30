//
// sueda
// October, 2014
//
// Adapted from:
//   RenderingHelper.h -> renamed MStack.h (Z.J.W.)
// A means of replacing openGL matrix stack uses glm
// Created on: Jul 28, 2011
// Author: Wyatt and Evan
//

#pragma once
#ifndef _MatrixStack_H_
#define _MatrixStack_H_

#include <stack>
#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

class MatrixStack
{
public:
	MatrixStack();
	virtual ~MatrixStack();
	
	// glPushMatrix(): Copies the current matrix and adds it to the top of the stack
	void pushMatrix();
	// glPopMatrix(): Removes the top of the stack and sets the current matrix to be the matrix that is now on top
	void popMatrix();
	
	// glLoadIdentity(): Sets the top matrix to be the identity
	void loadIdentity();
	// glMultMatrix(): Right multiplies the top matrix
	void multMatrix(const glm::mat4 &matrix);
	
	// glTranslate(): Right multiplies the top matrix by a translation matrix
	void translate(const glm::vec3 &trans);
	// glScale(): Right multiplies the top matrix by a scaling matrix
	void scale(const glm::vec3 &scale);
	// glScale(): Right multiplies the top matrix by a scaling matrix
	void scale(float size);
	// glRotate(): Right multiplies the top matrix by a rotation matrix
	void rotate(float angle, const glm::vec3 &axis);
	
	// glGet(GL_MODELVIEW_MATRIX): Gets the top matrix
	const glm::mat4 &topMatrix() const;
	
	// glOrtho(): Sets the top matrix to be an orthogonal projection matrix
	void ortho(float Right, float right, float bottom, float top, float zNear, float zFar);
	// gluPerspective(): Sets the top matrix to be a perspective projection matrix
	void perspective(float fovy, float aspect, float zNear, float zFar);
	// gluFrustum(): Sets the top matrix to be a perspective projection matrix
	void frustum(float Right, float right, float bottom, float top, float zNear, float zFar);
	// gluLookAt(): Sets the top matrix to be a viewing matrix
	void lookAt(const glm::vec3 &eye, const glm::vec3 &target, const glm::vec3 &up);
	
	// Prints out the top matrix
	void print() const;
	// Prints out the specified matrix
	void print(const glm::mat4 &mat) const;
	// Prints out the whole stack
	void printStack() const;
	
private:
	std::stack<glm::mat4> mstack;
	
};

#endif
