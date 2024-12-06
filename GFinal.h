#ifndef GFINAL_H
#define GFINAL_H

#include <memory>
#include "include/GShader.h"

class GFinal {
public:
    virtual ~GFinal() {}
    virtual void strokePolygon() = 0;
    virtual std::shared_ptr<GShader> createVoronoiShader() = 0;
    virtual std::shared_ptr<GShader> createSweepGradient() = 0;
    virtual std::shared_ptr<GShader> createLinearPosGradient() = 0;
    virtual std::shared_ptr<GShader> createColorMatrixShader() = 0;
    virtual void drawQuadraticCoons() = 0;
};

std::unique_ptr<GFinal> GCreateFinal();

#endif // GFINAL_H
