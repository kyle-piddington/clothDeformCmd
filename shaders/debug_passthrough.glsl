attribute vec4 vertPos; // in object space
uniform mat4 P;
uniform mat4 MV;
void main()
{
    gl_Position = P * MV * vertPos;
}