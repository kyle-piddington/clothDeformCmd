attribute vec4 vertPos; // in object space
attribute vec3 vertNor; // in object space
//attribute vec2 vertTexCoords;
uniform mat4 P;
uniform mat4 MV;
uniform mat4 NORM;

varying vec3 normal;
varying vec3 position;
void main()
{
    position = (MV*vertPos).xyz;
    normal = (NORM*vec4(vertNor,0.0)).xyz;
    gl_Position = P * vec4(position,1.0);
}
