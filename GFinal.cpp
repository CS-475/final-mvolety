#include "GFinal.h"
#include "include/GColor.h"

// Return null for shaders we don't implement
std::unique_ptr<GShader> MyFinal::createVoronoiShader() {
    return nullptr;
}

std::unique_ptr<GShader> MyFinal::createSweepGradient() {
    return nullptr;
}

// Return a solid red shader for demonstration
std::unique_ptr<GShader> MyFinal::createLinearPosGradient() {
    GColor red = {1, 1, 0, 0}; // Opaque red
    return std::make_unique<SolidShader>(red);
}

// Return a solid blue shader for demonstration
std::unique_ptr<GShader> MyFinal::createColorMatrixShader() {
    GColor blue = {1, 0, 0, 1}; // Opaque blue
    return std::make_unique<SolidShader>(blue);
}

void MyFinal::strokePolygon() {
    // no-op
}

void MyFinal::drawQuadraticCoons() {
    // no-op
}

std::unique_ptr<GFinal> GCreateFinal() {
    return std::make_unique<MyFinal>();
}
