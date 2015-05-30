//
// sueda
// October, 2014
//

#pragma  once
#ifndef __Camera__
#define __Camera__

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"

class MatrixStack;

class Camera
{
public:
	
	enum {
		ROTATE = 0,
		TRANSLATE,
		SCALE
	};
	
	Camera();
	virtual ~Camera();
	void setAspect(float a) { aspect = a; };
	void setRotationFactor(float f) { rfactor = f; };
	void setTranslationFactor(float f) { tfactor = f; };
	void setScaleFactor(float f) { sfactor = f; };
	void mouseClicked(int x, int y, bool shift, bool ctrl, bool alt);
	void mouseMoved(int x, int y);
	void applyProjectionMatrix(MatrixStack *P) const;
	void applyOrthoganalMatrix(MatrixStack *P, float xWdth) const;
	
	void applyCameraMatrix(MatrixStack *MV) const;
	
private:
	float aspect;
	float fovy;
	float znear;
	float zfar;
	glm::vec2 rotations;
	glm::vec3 translations;
	float scale;
	glm::vec2 mousePrev;
	int state;
	float rfactor;
	float tfactor;
	float sfactor;
};

#endif
