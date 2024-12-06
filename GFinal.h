#ifndef GFINAL_H
#define GFINAL_H

#include <memory>
#include "include/GShader.h"
#include "include/GColor.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"

// The GFinal base class
class GFinal {
public:
    virtual ~GFinal() {}

    virtual std::unique_ptr<GShader> createVoronoiShader() = 0;
    virtual std::unique_ptr<GShader> createSweepGradient() = 0;
    virtual std::unique_ptr<GShader> createLinearPosGradient() = 0;
    virtual std::unique_ptr<GShader> createColorMatrixShader() = 0;
    virtual void strokePolygon() = 0;
    virtual void drawQuadraticCoons() = 0;
};

// A very simple solid-color shader to avoid complexity
class SolidShader : public GShader {
public:
    SolidShader(GColor c) {
        GColor pinned = c.pinToUnit();
        float a = pinned.a, r = pinned.r * a, g = pinned.g * a, b = pinned.b * a;
        fPixel = GPixel_PackARGB((int)(a*255 + 0.5f),
                                 (int)(r*255 + 0.5f),
                                 (int)(g*255 + 0.5f),
                                 (int)(b*255 + 0.5f));
    }

    bool isOpaque() override {
        return GPixel_GetA(fPixel) == 0xFF;
    }
    bool setContext(const GMatrix&) override {
        return true;
    }
    void shadeRow(int x, int y, int count, GPixel row[]) override {
        for (int i = 0; i < count; ++i) {
            row[i] = fPixel;
        }
    }

private:
    GPixel fPixel;
};

// Our subclass that implements GFinal
class MyFinal : public GFinal {
public:
    std::unique_ptr<GShader> createVoronoiShader() override;
    std::unique_ptr<GShader> createSweepGradient() override;
    std::unique_ptr<GShader> createLinearPosGradient() override;
    std::unique_ptr<GShader> createColorMatrixShader() override;
    void strokePolygon() override;
    void drawQuadraticCoons() override;
};

std::unique_ptr<GFinal> GCreateFinal();

#endif // GFINAL_H
