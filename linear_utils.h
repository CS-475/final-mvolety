#ifndef LINEAR_UTILS_H
#define LINEAR_UTILS_H

#include "include/GShader.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include "include/GMatrix.h"

// The class definition goes here
class LinearGradientShader : public GShader {
public:
    LinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode);
    
    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    GColor interpolateColor(float t) const;
    GPixel colorToPixel(const GColor& color) const;
    float applyTileMode(float t) const;

    GPoint fP0, fP1;
    std::vector<GColor> fColors;
    GTileMode fTileMode;
    GMatrix fLocalMatrix;
};

#endif // LINEAR_UTILS_H
