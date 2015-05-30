//
// sueda
// October, 2014
//

#include "Camera.h"
#include "MatrixStack.h"
#include <iostream>

Camera::Camera() :
	aspect(1.0f),
	fovy(45.0f/180.0f*M_PI),
	znear(0.1f),
	zfar(1000.0f),
	rotations(0.0, 0.0),
	translations(0.0f, 0.0f, -2.0f),
	rfactor(0.01f),
	tfactor(0.005f),
	sfactor(0.005f)
{
}

Camera::~Camera()
{
	
}

void Camera::mouseClicked(int x, int y, bool shift, bool ctrl, bool alt)
{
	mousePrev.x = x;
	mousePrev.y = y;
	if(shift) {
		state = Camera::SCALE;
	} else if(ctrl) {
		state = Camera::TRANSLATE;
	} else {
		state = Camera::ROTATE;
	}
}

void Camera::mouseMoved(int x, int y)
{
	glm::vec2 mouseCurr(x, y);
	glm::vec2 dv = mouseCurr - mousePrev;
	switch(state) {
		case Camera::ROTATE:
			rotations += rfactor * dv;
			break;
		case Camera::TRANSLATE:
			translations.x += tfactor * dv.x;
			translations.y -= tfactor * dv.y;
			break;
		case Camera::SCALE:
			translations.z *= (1.0f + sfactor * dv.y);
			break;
	}
	mousePrev.x = x;
	mousePrev.y = y;
}

void Camera::applyProjectionMatrix(MatrixStack *P) const
{
	P->perspective(fovy, aspect, znear, zfar);
}

void Camera::applyOrthoganalMatrix(MatrixStack *P, float xWdth) const
{
	P->ortho(-xWdth/2.0,xWdth/2.0,(-xWdth/2.0)*aspect,(xWdth/2.0)*aspect,znear,zfar);
}
void Camera::applyCameraMatrix(MatrixStack *MV) const
{
	MV->translate(translations);
	MV->rotate(rotations.y, glm::vec3(1.0f, 0.0f, 0.0f));
	MV->rotate(rotations.x, glm::vec3(0.0f, 1.0f, 0.0f));
}
