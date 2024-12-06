#include "include/GPath.h"
#include "include/GMatrix.h"
#include <algorithm> 

GRect GPath::bounds() const {
    if (fPts.empty()) return GRect::XYWH(0, 0, 0, 0);

    float minX = fPts[0].x, minY = fPts[0].y;
    float maxX = fPts[0].x, maxY = fPts[0].y;

    for (const auto& pt : fPts) {
        minX = std::min(minX, pt.x);
        minY = std::min(minY, pt.y);
        maxX = std::max(maxX, pt.x);
        maxY = std::max(maxY, pt.y);
    }

    const float epsilon = 0.5f;
    return GRect::LTRB(minX - epsilon, minY - epsilon, maxX + epsilon, maxY + epsilon);
}



GPoint lerp(const GPoint& a, const GPoint& b, float t) {
    return { a.x + t * (b.x - a.x), a.y + t * (b.y - a.y) };
}

void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    GPoint ab = lerp(src[0], src[1], t);
    GPoint bc = lerp(src[1], src[2], t);
    GPoint abc = lerp(ab, bc, t);
    dst[0] = src[0];
    dst[1] = ab;
    dst[2] = abc;
    dst[3] = bc;
    dst[4] = src[2];
}

void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    GPoint ab = lerp(src[0], src[1], t);
    GPoint bc = lerp(src[1], src[2], t);
    GPoint cd = lerp(src[2], src[3], t);
    GPoint abc = lerp(ab, bc, t);
    GPoint bcd = lerp(bc, cd, t);
    GPoint abcd = lerp(abc, bcd, t);
    dst[0] = src[0];
    dst[1] = ab;
    dst[2] = abc;
    dst[3] = abcd;
    dst[4] = bcd;
    dst[5] = cd;
    dst[6] = src[3];
}
