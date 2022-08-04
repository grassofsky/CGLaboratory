const char *fragmentShaderSource = R"(

#version 420

layout(binding=0, offset=0) uniform atomic_uint Counter;
layout(binding=1, r32ui) uniform uimage2D HeadPointerImage;

out vec4 oFragColor;

void main()
{
    uint index = atomicCounterIncrement(Counter);
    uint oldHead = imageAtomicExchange(HeadPointerImage, ivec2(gl_FragCoord.xy), uint(index));
    
    oFragColor = vec4(float(index)/50.0, float(index)/50.0, float(index)/50.0, 1.0);
}
")