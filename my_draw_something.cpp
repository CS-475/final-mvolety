#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "include/GShader.h"
#include <cmath>
#include "shader_utils.h"
#include "my_utils.h"

std::string GDrawSomething(GCanvas* canvas, GISize dimension) {
     canvas->clear(GColor::RGBA(0.9f, 0.9f, 0.9f, 1.0f));

    float centerX = dimension.width / 2.0f;
    float centerY = dimension.height / 2.0f;

    GBitmap bitmap;
    GMatrix localMatrix; 

    const int circlePoints = 100;
    GPoint circle[circlePoints];
    float circleRadius = 100.0f;

    for (int i = 0; i < circlePoints; ++i) {
        float angle = 2.0f * static_cast<float>(M_PI) * i / circlePoints;
        circle[i] = {
            centerX + circleRadius * std::cos(angle),
            centerY + circleRadius * std::sin(angle)
        };
    }

    GPaint circlePaint(GColor::RGBA(0.3f, 0.5f, 0.9f, 1.0f));
    canvas->drawConvexPolygon(circle, circlePoints, circlePaint);

    return "painting";
}
