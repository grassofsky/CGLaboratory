#include <string>

#ifndef HEAD_POINTER_FRAG_H_
#define HEAD_POINTER_FRAG_H_

std::string fragmentShaderSource = R"(

#version 420

layout(binding=0, offset=0) uniform atomic_uint Counter;
layout(binding=1, r32ui) uniform uimage2D HeadPointerImage;

out vec4 oFragColor;

void main()
{
    uint index = atomicCounterIncrement(Counter);
    uint oldHead = imageAtomicExchange(HeadPointerImage, ivec2(gl_FragCoord.xy), uint(index));
    
    oFragColor = vec4(float(index)/480000.0, float(index)/480000.0, float(index)/480000., 1.0);
}
)";

#endif