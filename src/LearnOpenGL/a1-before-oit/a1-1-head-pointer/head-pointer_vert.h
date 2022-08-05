#include <string>

#ifndef HEAD_POINTER_VERT_H_
#define HEAD_POINTER_VERT_H_

std::string vertexShaderSource = R"(
#version 420

layout (location = 0) in vec3 vPos;
layout (location = 2) in vec3 vNormal;
layout (location = 3) in vec4 vColor;

uniform mat4 osg_ModelViewProjectionMatrix;

out vec3 v2fPos;
out vec3 v2fNormal;
out vec4 v2fColor;

void main()
{
    gl_Position = osg_ModelViewProjectionMatrix * vec4(vPos, 1.0);
}
)";

#endif
