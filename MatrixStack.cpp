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

#include "MatrixStack.h"
#include <stdio.h>
#include <vector>
#include "glm/gtx/transform.hpp"

MatrixStack::MatrixStack()
{
	mstack.push(glm::mat4(1.0f));
}

MatrixStack::~MatrixStack()
{
}

void MatrixStack::pushMatrix()
{
	const glm::mat4 &top = mstack.top();
	mstack.push(top);
	assert(mstack.size() < 100);
}

void MatrixStack::popMatrix()
{
	assert(!mstack.empty());
	mstack.pop();
	// There should always be one matrix left.
	assert(!mstack.empty());
}

void MatrixStack::loadIdentity()
{
	glm::mat4 &top = mstack.top();
	top = glm::mat4(1.0f);
}

void MatrixStack::translate(const glm::vec3 &trans)
{
	glm::mat4 &top = mstack.top();
	top *= glm::translate(trans);
}

void MatrixStack::scale(const glm::vec3 &scale)
{
	glm::mat4 &top = mstack.top();
	top *= glm::scale(scale);
}

void MatrixStack::scale(float s)
{
	glm::mat4 &top = mstack.top();
	top *= glm::scale(glm::vec3(s, s, s));
}

void MatrixStack::rotate(float angle, const glm::vec3 &axis)
{
	glm::mat4 &top = mstack.top();
	top *= glm::rotate(angle, axis);
}

void MatrixStack::multMatrix(const glm::mat4 &matrix)
{
	glm::mat4 &top = mstack.top();
	top *= matrix;
}

void MatrixStack::ortho(float l, float r, float b, float t, float zNear, float zFar)
{
	glm::mat4 &top = mstack.top();
	top = glm::ortho(l, r, b, t, zNear, zFar);
}

void MatrixStack::perspective(float fovy, float aspect, float zNear, float zFar)
{
	glm::mat4 &top = mstack.top();
	top = glm::perspective(fovy, aspect, zNear, zFar);
}

void MatrixStack::frustum(float l, float r, float b, float t, float zNear, float zFar)
{
	glm::mat4 &top = mstack.top();
	top = glm::frustum(l, r, b, t, zNear, zFar);
}

void MatrixStack::lookAt(const glm::vec3 &eye, const glm::vec3 &target, const glm::vec3 &up)
{
	glm::mat4 &top = mstack.top();
	top = glm::lookAt(eye, target, up);
}

const glm::mat4 &MatrixStack::topMatrix() const
{
	return mstack.top();
}

void MatrixStack::print(const glm::mat4 &mat) const
{
	// Prints so that translation is on the last column
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%- 5.2f ", mat[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void MatrixStack::print() const
{
	// Prints so that translation is on the last column
	const glm::mat4 &mat = mstack.top();
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%- 5.2f ", mat[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void MatrixStack::printStack() const
{
	// Copy everything to a non-const stack
	std::stack<glm::mat4> tempStack = mstack;
	while(!tempStack.empty()) {
		glm::mat4 &top = tempStack.top();
		print(top);
		tempStack.pop();
	}
}
