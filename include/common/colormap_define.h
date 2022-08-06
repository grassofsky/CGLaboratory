#ifndef COLORMAP_DEINFE_H_
#define COLORMAP_DEFINE_H_

#include <vector>

void InterplateColor(double value, double &red, double &green, double &blue)
{
    static std::vector<unsigned char> colormap = {
        12, 12, 242,
        12, 12, 242,
        12, 16, 242,
        12, 20, 242,
        12, 24, 242,
        12, 28, 242,
        12, 32, 242,
        12, 36, 242,
        12, 40, 242,
        12, 44, 242,
        12, 48, 242,
        12, 52, 242,
        12, 56, 242,
        12, 59, 242,
        12, 63, 242,
        12, 67, 242,
        12, 70, 242,
        12, 74, 242,
        12, 78, 242,
        12, 82, 242,
        12, 86, 242,
        12, 90, 242,
        12, 94, 242,
        12, 98, 242,
        12, 102, 242,
        12, 106, 242,
        12, 110, 242,
        12, 114, 242,
        12, 118, 242,
        12, 122, 242,
        12, 126, 242,
        12, 130, 242,
        12, 133, 242,
        12, 137, 242,
        12, 140, 242,
        12, 143, 242,
        12, 147, 242,
        12, 150, 242,
        12, 154, 242,
        12, 157, 242,
        12, 161, 242,
        12, 164, 242,
        12, 168, 242,
        12, 172, 242,
        12, 175, 242,
        12, 179, 242,
        12, 182, 242,
        12, 186, 242,
        12, 189, 242,
        12, 193, 242,
        12, 196, 242,
        12, 200, 242,
        12, 203, 242,
        12, 207, 242,
        12, 210, 242,
        12, 214, 242,
        12, 217, 242,
        12, 221, 242,
        12, 224, 242,
        12, 228, 242,
        12, 231, 242,
        12, 235, 242,
        12, 238, 242,
        12, 242, 238,
        12, 242, 235,
        12, 242, 231,
        12, 242, 228,
        12, 242, 224,
        12, 242, 221,
        12, 242, 217,
        12, 242, 214,
        12, 242, 210,
        12, 242, 207,
        12, 242, 203,
        12, 242, 200,
        12, 242, 196,
        12, 242, 193,
        12, 242, 189,
        12, 242, 186,
        12, 242, 182,
        12, 242, 179,
        12, 242, 175,
        12, 242, 172,
        12, 242, 168,
        12, 242, 164,
        12, 242, 161,
        12, 242, 157,
        12, 242, 154,
        12, 242, 150,
        12, 242, 147,
        12, 242, 143,
        12, 242, 140,
        12, 242, 137,
        12, 242, 133,
        12, 242, 130,
        12, 242, 127,
        12, 242, 124,
        12, 242, 120,
        12, 242, 116,
        12, 242, 112,
        12, 242, 108,
        12, 242, 104,
        12, 242, 100,
        12, 242, 96,
        12, 242, 92,
        12, 242, 88,
        12, 242, 85,
        12, 242, 81,
        12, 242, 77,
        12, 242, 74,
        12, 242, 70,
        12, 242, 67,
        12, 242, 63,
        12, 242, 59,
        12, 242, 56,
        12, 242, 52,
        12, 242, 49,
        12, 242, 45,
        12, 242, 42,
        12, 242, 38,
        12, 242, 35,
        12, 242, 31,
        12, 242, 28,
        12, 242, 24,
        12, 242, 21,
        12, 242, 18,
        12, 242, 15,
        12, 242, 12,
        15, 242, 12,
        18, 242, 12,
        21, 242, 12,
        24, 242, 12,
        28, 242, 12,
        31, 242, 12,
        35, 242, 12,
        38, 242, 12,
        42, 242, 12,
        45, 242, 12,
        49, 242, 12,
        52, 242, 12,
        56, 242, 12,
        59, 242, 12,
        63, 242, 12,
        67, 242, 12,
        70, 242, 12,
        74, 242, 12,
        77, 242, 12,
        81, 242, 12,
        85, 242, 12,
        88, 242, 12,
        92, 242, 12,
        96, 242, 12,
        100, 242, 12,
        104, 242, 12,
        108, 242, 12,
        112, 242, 12,
        116, 242, 12,
        120, 242, 12,
        124, 242, 12,
        127, 242, 12,
        130, 242, 12,
        133, 242, 12,
        137, 242, 12,
        140, 242, 12,
        143, 242, 12,
        147, 242, 12,
        150, 242, 12,
        154, 242, 12,
        157, 242, 12,
        161, 242, 12,
        164, 242, 12,
        168, 242, 12,
        172, 242, 12,
        175, 242, 12,
        179, 242, 12,
        182, 242, 12,
        186, 242, 12,
        189, 242, 12,
        193, 242, 12,
        196, 242, 12,
        200, 242, 12,
        203, 242, 12,
        207, 242, 12,
        210, 242, 12,
        214, 242, 12,
        217, 242, 12,
        221, 242, 12,
        224, 242, 12,
        228, 242, 12,
        231, 242, 12,
        235, 242, 12,
        238, 242, 12,
        242, 242, 12,
        242, 238, 12,
        242, 235, 12,
        242, 231, 12,
        242, 228, 12,
        242, 224, 12,
        242, 221, 12,
        242, 217, 12,
        242, 214, 12,
        242, 210, 12,
        242, 207, 12,
        242, 203, 12,
        242, 200, 12,
        242, 196, 12,
        242, 193, 12,
        242, 189, 12,
        242, 186, 12,
        242, 182, 12,
        242, 179, 12,
        242, 175, 12,
        242, 172, 12,
        242, 168, 12,
        242, 164, 12,
        242, 161, 12,
        242, 157, 12,
        242, 154, 12,
        242, 150, 12,
        242, 147, 12,
        242, 143, 12,
        242, 140, 12,
        242, 137, 12,
        242, 133, 12,
        242, 130, 12,
        242, 126, 12,
        242, 122, 12,
        242, 118, 12,
        242, 114, 12,
        242, 110, 12,
        242, 106, 12,
        242, 102, 12,
        242, 98, 12,
        242, 94, 12,
        242, 90, 12,
        242, 86, 12,
        242, 82, 12,
        242, 78, 12,
        242, 74, 12,
        242, 70, 12,
        242, 67, 12,
        242, 63, 12,
        242, 59, 12,
        242, 56, 12,
        242, 52, 12,
        242, 48, 12,
        242, 44, 12,
        242, 40, 12,
        242, 36, 12,
        242, 32, 12,
        242, 28, 12,
        242, 24, 12,
        242, 20, 12,
        242, 16, 12,
        242, 12, 12,
        242, 12, 12};

    int colorSize = colormap.size() / 3;
    if (value > 1.0) value = 1.0;
    if (value < 0) value = 0;

    double dIdx = value * (colorSize - 1);
    int idxBefore = static_cast<int>(dIdx);
    int idxAfter = idxBefore + 1;
    if (idxAfter > colorSize-1)
    {
        idxAfter = colorSize - 1;
    }

    double factor = dIdx - idxBefore;
    unsigned char before[] = {colormap[idxBefore*3], colormap[idxBefore*3+1], colormap[idxBefore*3+2]};
    unsigned char after[] = {colormap[idxAfter*3], colormap[idxAfter*3+1], colormap[idxAfter*3+2]};
    
    red = ((1-factor) * before[0] + factor * after[0])/255.0;
    green = ((1-factor) * before[1] + factor * after[1])/255.0;
    blue = ((1-factor) * before[2] + factor * after[2])/255.0;
}

#endif