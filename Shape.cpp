//
// sueda
// October, 2014
//

#include <iostream>
#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "Shape.h"
#include "GLSL.h"

using namespace std;

Shape::Shape() :
	posBufID(0),
	norBufID(0),
	indBufID(0)
{
}

Shape::~Shape()
{
}

void Shape::load(const string &meshName)
{
	// Load geometry
	// Some obj files contain material information.
	// We'll ignore them for this assignment.
	std::vector<tinyobj::material_t> objMaterials;
	string err = tinyobj::LoadObj(shapes, objMaterials, meshName.c_str());
	if(!err.empty()) {
		cerr << err << endl;
	}
	
	// Scale the vertex positions so that they fit within [-1, +1] in all three dimensions.
	vector<float> &posBuf = shapes[0].mesh.positions;
	glm::vec3 vmin(posBuf[0], posBuf[1], posBuf[2]);
	glm::vec3 vmax(posBuf[0], posBuf[1], posBuf[2]);
	for(int i = 3; i < (int)posBuf.size(); i += 3) {
		glm::vec3 v(posBuf[i], posBuf[i+1], posBuf[i+2]);
		vmin.x = min(vmin.x, v.x);
		vmin.y = min(vmin.y, v.y);
		vmin.z = min(vmin.z, v.z);
		vmax.x = max(vmax.x, v.x);
		vmax.y = max(vmax.y, v.y);
		vmax.z = max(vmax.z, v.z);
	}
	glm::vec3 center = 0.5f*(vmin + vmax);
	glm::vec3 diff = vmax - vmin;
	float diffmax = diff.x;
	diffmax = max(diffmax, diff.y);
	diffmax = max(diffmax, diff.z);
	float scale = 1.0f / diffmax;
	for(int i = 0; i < (int)posBuf.size(); i += 3) {
		posBuf[i  ] = (posBuf[i  ] - center.x) * scale;
		posBuf[i+1] = (posBuf[i+1] - center.y) * scale;
		posBuf[i+2] = (posBuf[i+2] - center.z) * scale;
	}
}

void Shape::init()
{
	// Send the position array to the GPU
	const vector<float> &posBuf = shapes[0].mesh.positions;
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STATIC_DRAW);
	
	// Send the normal array (if it exists) to the GPU
	const vector<float> &norBuf = shapes[0].mesh.normals;
	if(!norBuf.empty()) {
		glGenBuffers(1, &norBufID);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_STATIC_DRAW);
	} else {
		norBufID = 0;
	}
	
	// Send the index array to the GPU
	const vector<unsigned int> &indBuf = shapes[0].mesh.indices;
	glGenBuffers(1, &indBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indBuf.size()*sizeof(unsigned int), &indBuf[0], GL_STATIC_DRAW);
	
	// Unbind the arrays
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	assert(glGetError() == GL_NO_ERROR);
}

void Shape::draw(GLint h_pos, GLint h_nor)
{
	// Enable and bind position array for drawing
	GLSL::enableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, 0);
	
	// Enable and bind normal array (if it exists) for drawing
	if(norBufID) {
		GLSL::enableVertexAttribArray(h_nor);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, 0);
	}
	
	// Bind index array for drawing
	int nIndices = (int)shapes[0].mesh.indices.size();
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
	
	// Draw
	glDrawElements(GL_TRIANGLES, nIndices, GL_UNSIGNED_INT, 0);
	
	// Disable and unbind
	if(norBufID) {
		GLSL::disableVertexAttribArray(h_nor);
	}
	GLSL::disableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
