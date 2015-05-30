uniform vec3 lightPos[10];
uniform float lightIntensity[10];
uniform int numLights;
uniform vec3 ka;
uniform vec3 kd;
uniform vec3 ks;
uniform float s;
uniform mat3 T;
uniform sampler2D texture;

varying vec3 normal; // passed from the vertex shader
varying vec3 position;
varying vec2 fragTexCoords;

void main()
{
    vec3 color = vec3(0);
    vec3 n = normalize(normal);
    for(int i = 0; i < numLights; i++){
        vec3 h = normalize(
            normalize(lightPos[i]-position)-normalize(position));
          color+= lightIntensity[i]*
          (ka + kd*max(dot(normalize(lightPos[i]-position),n),0.0) +
            ks*pow(max(dot(h,n),0.0),s));
	}
    gl_FragColor = vec4(color.r, color.g, color.b, 1.0);
}
