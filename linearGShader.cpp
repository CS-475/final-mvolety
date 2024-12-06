#include "linear_utils.h"
#include <cmath>
#include <algorithm>

LinearGradientShader::LinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode)
    : fP0(p0), fP1(p1), fColors(colors, colors + count), fTileMode(mode) {}

bool LinearGradientShader::isOpaque() {
    return std::all_of(fColors.begin(), fColors.end(), [](const GColor& color) { return color.a >= 1.0f; });
}

bool LinearGradientShader::setContext(const GMatrix& ctm) {
    float dx = fP1.x - fP0.x;
    float dy = fP1.y - fP0.y;
    float D = std::sqrt(dx * dx + dy * dy);
    if (D == 0) return false;

    GMatrix gradientMatrix(dx, -dy, 0, dy, dx, 0);
    GMatrix translationMatrix(1, 0, fP0.x, 0, 1, fP0.y);
    GMatrix combinedMatrix = GMatrix::Concat(translationMatrix, gradientMatrix);
    combinedMatrix = GMatrix::Concat(ctm, combinedMatrix);

    auto inv = combinedMatrix.invert();
    if (!inv.has_value()) return false;
    fLocalMatrix = inv.value();
    return true;
}

void LinearGradientShader::shadeRow(int x, int y, int count, GPixel row[]) {
    for (int i = 0; i < count; ++i) {
        GPoint devicePoint = {x + i + 0.5f, y + 0.5f};
        GPoint localPoint;
        fLocalMatrix.mapPoints(&localPoint, &devicePoint, 1);

        float t = applyTileMode(localPoint.x);
        GColor color = interpolateColor(t);
        row[i] = colorToPixel(color);
    }
}

GColor LinearGradientShader::interpolateColor(float t) const {
    int n = fColors.size();
    int idx = static_cast<int>(t * (n - 1));
    float localT = (t * (n - 1)) - idx;
    GColor c0 = fColors[std::min(idx, n - 1)];
    GColor c1 = fColors[std::min(idx + 1, n - 1)];
    return {
        c0.r * (1 - localT) + c1.r * localT,
        c0.g * (1 - localT) + c1.g * localT,
        c0.b * (1 - localT) + c1.b * localT,
        c0.a * (1 - localT) + c1.a * localT
    };
}

GPixel LinearGradientShader::colorToPixel(const GColor& color) const {
    uint8_t a = static_cast<uint8_t>(color.a * 255 + 0.5f);
    uint8_t r = static_cast<uint8_t>(color.r * color.a * 255 + 0.5f);
    uint8_t g = static_cast<uint8_t>(color.g * color.a * 255 + 0.5f);
    uint8_t b = static_cast<uint8_t>(color.b * color.a * 255 + 0.5f);
    return GPixel_PackARGB(a, r, g, b);
}

float LinearGradientShader::applyTileMode(float t) const {
    switch (fTileMode) {
        case GTileMode::kClamp:
            return std::clamp(t, 0.0f, 1.0f);
        case GTileMode::kRepeat:
            return t - std::floor(t);  // Keep t within the [0, 1) range
        case GTileMode::kMirror: {
            float modT = std::fabs(t - std::floor(t));  // Mirror t between [0, 1)
            return (static_cast<int>(std::floor(t)) % 2 == 0) ? modT : 1.0f - modT;
        }
        default:
            return t;
    }
}

std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode) {
    if (count < 2) {
        return nullptr;  
    }
    return std::make_shared<LinearGradientShader>(p0, p1, colors, count, mode);
}
