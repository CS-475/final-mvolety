#include "my_utils.h"
#include "include/GShader.h" 
#include "include/GPathBuilder.h"
#include "shader_utils.h"
#include "include/GPath.h"
#include <algorithm>
#include <cmath>

GPixel convertColorToPix(const GColor& color) {
    GColor clampedColor = color.pinToUnit();
    uint8_t a = static_cast<uint8_t>(clampedColor.a * 255 + 0.5f);
    uint8_t r = static_cast<uint8_t>(clampedColor.r * clampedColor.a * 255 + 0.5f);
    uint8_t g = static_cast<uint8_t>(clampedColor.g * clampedColor.a * 255 + 0.5f);
    uint8_t b = static_cast<uint8_t>(clampedColor.b * clampedColor.a * 255 + 0.5f);
    return GPixel_PackARGB(a, r, g, b);
}

GColor interpolateColor(const GColor& c0, const GColor& c1, const GColor& c2, const GColor& c3,
                        float w0, float w1, float w2, float w3) {
    GColor result;
    result.a = std::clamp(w0 * c0.a + w1 * c1.a + w2 * c2.a + w3 * c3.a, 0.0f, 1.0f);
    result.r = std::clamp(w0 * c0.r + w1 * c1.r + w2 * c2.r + w3 * c3.r, 0.0f, 1.0f);
    result.g = std::clamp(w0 * c0.g + w1 * c1.g + w2 * c2.g + w3 * c3.g, 0.0f, 1.0f);
    result.b = std::clamp(w0 * c0.b + w1 * c1.b + w2 * c2.b + w3 * c3.b, 0.0f, 1.0f);
    return result;
}

void MyCanvas::clear(const GColor& color) {
    GPixel pixColor = convertColorToPix(color);
    int height = fDevice.height();
    int width = fDevice.width();

    for (int y = 0; y < height; y++) {
        GPixel* rowAddr = fDevice.getAddr(0, y);
        std::fill(rowAddr, rowAddr + width, pixColor);
    }
}

void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
                        int count, const int indices[], const GPaint& paint) {
    for (int i = 0; i < count; ++i) {
        int idx0 = indices[i * 3];
        int idx1 = indices[i * 3 + 1];
        int idx2 = indices[i * 3 + 2];

        const GColor* color0 = colors ? &colors[idx0] : nullptr;
        const GColor* color1 = colors ? &colors[idx1] : nullptr;
        const GColor* color2 = colors ? &colors[idx2] : nullptr;

        const GPoint* tex0 = texs ? &texs[idx0] : nullptr;
        const GPoint* tex1 = texs ? &texs[idx1] : nullptr;
        const GPoint* tex2 = texs ? &texs[idx2] : nullptr;

        drawTriangle(verts[idx0], verts[idx1], verts[idx2],
                     color0, color1, color2,
                     tex0, tex1, tex2,
                     paint);
    }
}

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                        int level, const GPaint& paint) {
    int subdivisions = level + 1;
    int gridSize = subdivisions + 1;

    std::vector<GPoint> meshVerts(gridSize * gridSize);
    std::vector<GColor> meshColors;
    std::vector<GPoint> meshTexs;

    if (colors) meshColors.resize(meshVerts.size());
    if (texs) meshTexs.resize(meshVerts.size());

    for (int i = 0; i <= subdivisions; ++i) {
        float t = static_cast<float>(i) / subdivisions;
        for (int j = 0; j <= subdivisions; ++j) {
            float s = static_cast<float>(j) / subdivisions;
            int idx = i * gridSize + j;

            meshVerts[idx] = {
                verts[0].x * (1 - s) * (1 - t) + verts[1].x * s * (1 - t) +
                verts[2].x * s * t + verts[3].x * (1 - s) * t,
                verts[0].y * (1 - s) * (1 - t) + verts[1].y * s * (1 - t) +
                verts[2].y * s * t + verts[3].y * (1 - s) * t
            };

            if (colors) {
                meshColors[idx] = interpolateColor(
                    colors[0], colors[1], colors[2], colors[3],
                     (1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t
                );

            }

            if (texs) {
                meshTexs[idx] = {
                    texs[0].x * (1 - s) * (1 - t) + texs[1].x * s * (1 - t) +
                    texs[2].x * s * t + texs[3].x * (1 - s) * t,
                    texs[0].y * (1 - s) * (1 - t) + texs[1].y * s * (1 - t) +
                    texs[2].y * s * t + texs[3].y * (1 - s) * t
                };
            }
        }
    }

    std::vector<int> indices;
    for (int i = 0; i < subdivisions; ++i) {
        for (int j = 0; j < subdivisions; ++j) {
            int idx0 = i * gridSize + j;
            int idx1 = idx0 + 1;
            int idx2 = idx0 + gridSize + 1;
            int idx3 = idx0 + gridSize;

            indices.push_back(idx0);
            indices.push_back(idx1);
            indices.push_back(idx2);

            indices.push_back(idx0);
            indices.push_back(idx2);
            indices.push_back(idx3);
        }
    }

    drawMesh(meshVerts.data(),
             colors ? meshColors.data() : nullptr,
             texs ? meshTexs.data() : nullptr,
             indices.size() / 3,
             indices.data(),
             paint);
}


void MyCanvas::computeBarycentric(int x, int y, const GPoint& p0, const GPoint& p1, const GPoint& p2,
                                  float& alpha, float& beta, float& gamma) {
    float denom = (p1.y - p2.y) * (p0.x - p2.x) + (p2.x - p1.x) * (p0.y - p2.y);
    if (std::abs(denom) < 1e-6f) {
        alpha = beta = gamma = -1;  
        return;
    }

    float invDenom = 1.0f / denom;
    alpha = ((p1.y - p2.y) * (x - p2.x) + (p2.x - p1.x) * (y - p2.y)) * invDenom;
    beta = ((p2.y - p0.y) * (x - p2.x) + (p0.x - p2.x) * (y - p2.y)) * invDenom;
    gamma = 1.0f - alpha - beta;
}


void MyCanvas::drawTriangle(const GPoint& p0, const GPoint& p1, const GPoint& p2,
                            const GColor* c0, const GColor* c1, const GColor* c2,
                            const GPoint* t0, const GPoint* t1, const GPoint* t2,
                            const GPaint& paint) {

    GPoint pts[3] = {p0, p1, p2};
    fMatrix.mapPoints(pts, pts, 3);

    float minX = std::min({pts[0].x, pts[1].x, pts[2].x});
    float maxX = std::max({pts[0].x, pts[1].x, pts[2].x});
    float minY = std::min({pts[0].y, pts[1].y, pts[2].y});
    float maxY = std::max({pts[0].y, pts[1].y, pts[2].y});

    int x0 = std::max(0, static_cast<int>(std::floor(minX)));
    int x1 = std::min(fDevice.width() - 1, static_cast<int>(std::ceil(maxX)));
    int y0 = std::max(0, static_cast<int>(std::floor(minY)));
    int y1 = std::min(fDevice.height() - 1, static_cast<int>(std::ceil(maxY)));

    float edge0_x = pts[1].x - pts[0].x, edge0_y = pts[1].y - pts[0].y;
    float edge1_x = pts[2].x - pts[1].x, edge1_y = pts[2].y - pts[1].y;
    float edge2_x = pts[0].x - pts[2].x, edge2_y = pts[0].y - pts[2].y;

    float area = edge0_x * edge2_y - edge0_y * edge2_x;
    if (std::abs(area) < 1e-5) return; 
    float invArea = 1.0f / area;

    GColor defaultColor = paint.getColor();
    GShader* shader = paint.peekShader();
    bool useShader = shader && shader->setContext(fMatrix);
    BlendProc blendProc = findBlend(paint.getBlendMode());

    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            float px = x + 0.5f, py = y + 0.5f;

            float w0 = ((pts[1].x - pts[2].x) * (py - pts[2].y) - 
                        (pts[1].y - pts[2].y) * (px - pts[2].x)) * invArea;
            float w1 = ((pts[2].x - pts[0].x) * (py - pts[0].y) - 
                        (pts[2].y - pts[0].y) * (px - pts[0].x)) * invArea;
            float w2 = 1.0f - w0 - w1;

            const float EPSILON = -1e-5f;
            if (w0 >= EPSILON && w1 >= EPSILON && w2 >= EPSILON) {
                GPoint localCoord;
                if (t0 && t1 && t2) {
                    localCoord = {
                        w0 * t0->x + w1 * t1->x + w2 * t2->x,
                        w0 * t0->y + w1 * t1->y + w2 * t2->y
                    };
                } else {
                    localCoord = {px, py};
                }

                GColor finalColor = defaultColor;
                if (c0 && c1 && c2) {
                    finalColor = {
                        w0 * c0->r + w1 * c1->r + w2 * c2->r,
                        w0 * c0->g + w1 * c1->g + w2 * c2->g,
                        w0 * c0->b + w1 * c1->b + w2 * c2->b,
                        w0 * c0->a + w1 * c1->a + w2 * c2->a
                    };
                }

                GPixel srcPixel = convertColorToPix(finalColor);
                if (useShader) {
                    shader->shadeRow(localCoord.x, localCoord.y, 1, &srcPixel);
                }

                GPixel* dstPixel = fDevice.getAddr(x, y);
                *dstPixel = blendProc(srcPixel, *dstPixel);
            }
        }
    }
}


void MyCanvas::drawRect(const GRect& rect, const GPaint& paint) {
    if (rect.isEmpty()) {
        return;
    }

    GPoint points[4] = {
        {rect.left, rect.top},       // No parentheses
        {rect.right, rect.top},
        {rect.right, rect.bottom},
        {rect.left, rect.bottom}
    };

    GPoint transformedPoints[4];
    fMatrix.mapPoints(transformedPoints, points, 4);

    GIRect ir = GRect::LTRB(
        std::floor(std::min({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x, transformedPoints[3].x})),
        std::floor(std::min({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y, transformedPoints[3].y})),
        std::ceil(std::max({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x, transformedPoints[3].x})),
        std::ceil(std::max({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y, transformedPoints[3].y}))
    ).roundOut();

    int left = std::max(0, ir.left);   // No parentheses
    int top = std::max(0, ir.top);
    int right = std::min(fDevice.width(), ir.right);
    int bottom = std::min(fDevice.height(), ir.bottom);

    if (left >= right || top >= bottom) {
        return;
    }

    if (ir.width() <= 0 || ir.height() <= 0) {
    return;
    }

    GShader* shader = paint.peekShader();
    GPixel srcPixel;
    bool useShader = (shader != nullptr && shader->setContext(fMatrix));

    BlendProc proc = findBlend(paint.getBlendMode());

    for (int y = top; y < bottom; ++y) {
        GPixel* rowAddr = fDevice.getAddr(0, y);
        for (int x = left; x < right; ++x) {
            if (useShader) {
                shader->shadeRow(x, y, 1, &srcPixel);
            } else {
                srcPixel = convertColorToPix(paint.getColor());
            }
            rowAddr[x] = proc(srcPixel, rowAddr[x]);
        }
    }
}



GRect computePolygonBounds(const GPoint points[], int count) {
    if (count == 0) {
        return GRect::LTRB(0, 0, 0, 0);
    }

    float minX = points[0].x;
    float minY = points[0].y;
    float maxX = points[0].x;
    float maxY = points[0].y;

    for (int i = 1; i < count; ++i) {
        if (points[i].x < minX) minX = points[i].x;
        if (points[i].x > maxX) maxX = points[i].x;
        if (points[i].y < minY) minY = points[i].y;
        if (points[i].y > maxY) maxY = points[i].y;
    }

    return GRect::LTRB(minX, minY, maxX, maxY);
}


void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {
    if (count < 3) return;

    GPoint transformedPoints[count];
    fMatrix.mapPoints(transformedPoints, points, count);

    GRect bounds = computePolygonBounds(transformedPoints, count);
    if (bounds.isEmpty()) return;
    int minYInt = std::max(0, (int)std::ceil(bounds.top));
    int maxYInt = std::min(fDevice.height() - 1, (int)std::floor(bounds.bottom));

    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr && shader->setContext(fMatrix));
    BlendProc proc = findBlend(paint.getBlendMode());
    GPixel srcPixel;

    for (int y = minYInt; y <= maxYInt; ++y) {
        std::vector<float> intersects;
        for (int i = 0; i < count; ++i) {
            GPoint p0 = transformedPoints[i];
            GPoint p1 = transformedPoints[(i + 1) % count];
            if ((p0.y <= y && p1.y > y) || (p0.y > y && p1.y <= y)) {
                float t = (y - p0.y) / (p1.y - p0.y);
                intersects.push_back(p0.x + t * (p1.x - p0.x));
            }
        }
        std::sort(intersects.begin(), intersects.end());

        for (size_t i = 0; i < intersects.size(); i += 2) {
            int xMin = std::max(0, (int)std::ceil(intersects[i]));
            int xMax = std::min(fDevice.width() - 1, (int)std::floor(intersects[i + 1]));

            GPixel* row = fDevice.getAddr(0, y);
            for (int x = xMin; x <= xMax; ++x) {
                if (useShader) {
                    shader->shadeRow(x, y, 1, &srcPixel);
                } else {
                    srcPixel = convertColorToPix(paint.getColor());
                }
                row[x] = proc(srcPixel, row[x]);
            }
        }
    }
}

int computeQuadSegments(const GPoint pts[3], float tolerance) {
    float ax = pts[0].x - 2 * pts[1].x + pts[2].x;
    float ay = pts[0].y - 2 * pts[1].y + pts[2].y;
    float maxDist = std::sqrt(ax * ax + ay * ay);
    return static_cast<int>(std::ceil(std::sqrt(maxDist / tolerance)));

}

int computeCubicSegments(const GPoint pts[4], float tolerance) {
    float ax = -pts[0].x + 3 * (pts[1].x - pts[2].x) + pts[3].x;
    float ay = -pts[0].y + 3 * (pts[1].y - pts[2].y) + pts[3].y;
    float maxDist = std::sqrt(ax * ax + ay * ay);
    return static_cast<int>(std::ceil(std::sqrt(maxDist / tolerance)));
}

GPoint evalQuad(const GPoint pts[3], float t) {
    float u = 1 - t;
    return {
        u * u * pts[0].x + 2 * u * t * pts[1].x + t * t * pts[2].x,
        u * u * pts[0].y + 2 * u * t * pts[1].y + t * t * pts[2].y
    };
}

GPoint evalCubic(const GPoint pts[4], float t) {
    float u = 1 - t;
    return {
        u * u * u * pts[0].x + 3 * u * u * t * pts[1].x + 3 * u * t * t * pts[2].x + t * t * t * pts[3].x,
        u * u * u * pts[0].y + 3 * u * u * t * pts[1].y + 3 * u * t * t * pts[2].y + t * t * t * pts[3].y
    };
}

void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
    GPathBuilder builder;
    builder.reset();
    GPath::Iter iter(path);
    GPoint pts[GPath::kMaxNextPoints];
    const float tolerance = 0.25f;

    while (auto verb = iter.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kMove:
                builder.moveTo(pts[0]);
                break;
            case GPathVerb::kLine:
                builder.lineTo(pts[1]);
                break;
            case GPathVerb::kQuad: {
                int segments = computeQuadSegments(pts, tolerance);
                for (int i = 0; i < segments; ++i) {
                    float t1 = static_cast<float>(i + 1) / segments;
                    builder.lineTo(evalQuad(pts, t1));
                }
                break;
            }
            case GPathVerb::kCubic: {
                int segments = computeCubicSegments(pts, tolerance);
                for (int i = 0; i < segments; ++i) {
                    float t1 = static_cast<float>(i + 1) / segments;
                    builder.lineTo(evalCubic(pts, t1));
                }
                break;
            }
            default:
                break;
        }
    }

    builder.transform(fMatrix);
    auto transformedPath = builder.detach();
    
    GPath::Edger edger(*transformedPath);
    std::vector<Edge> edges;
    while (auto verb = edger.next(pts)) {
        if (verb.value() == GPathVerb::kLine) {
            Edge edge = makeEdge(pts[0], pts[1]);
            if (edge.isUseful()) {
                edges.push_back(edge);
            }
        }
    }
    
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.top < b.top || (a.top == b.top && a.currX < b.currX);
    });

    int canvasHeight = fDevice.height();
    int canvasWidth = fDevice.width();
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr && shader->setContext(fMatrix));
    GPixel srcPixel;

    for (int y = 0; y < canvasHeight; ++y) {
        std::vector<int> xIntervals;
        int winding = 0;

        for (auto& edge : edges) {
            if (edge.isValid(y)) {
                int x = static_cast<int>(std::round(edge.computeX(y)));
                if (x < 0 || x >= canvasWidth) continue;
                xIntervals.push_back(x);
                winding += edge.winding;
            }
        }

        std::sort(xIntervals.begin(), xIntervals.end());

        for (size_t i = 0; i + 1 < xIntervals.size(); i += 2) {
            int left = xIntervals[i];
            int right = xIntervals[i + 1];
            left = std::max(0, left);
            right = std::min(canvasWidth, right);

            for (int xSpan = left; xSpan < right; ++xSpan) {
                if (useShader) {
                    shader->shadeRow(xSpan, y, 1, &srcPixel);
                } else {
                    srcPixel = convertColorToPix(paint.getColor());
                }

                GPixel* dst = fDevice.getAddr(xSpan, y);
                *dst = blendPixel(*dst, srcPixel, paint.getBlendMode());
            }
        }
    }
}


void MyCanvas::save() {
    fMatrixStack.push(fMatrix);
}

void MyCanvas::restore() {
    if (!fMatrixStack.empty()) {
        fMatrix = fMatrixStack.top();
        fMatrixStack.pop();
    }
}

void MyCanvas::concat(const GMatrix& matrix) {
    fMatrix = GMatrix::Concat(fMatrix, matrix);
}

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& bitmap) {
    return std::unique_ptr<GCanvas>(new MyCanvas(bitmap));
}

using BlendFunc = GPixel (*)(GPixel src, GPixel dst);

typedef GPixel (*BlendProc)(GPixel, GPixel);

const BlendProc gProcs[] = {
    clear_mode,   
    src_mode,     
    dst_mode,     
    src_over_mode, 
    dst_over_mode,
    src_in_mode,   
    dst_in_mode,  
    src_out_mode, 
    dst_out_mode, 
    src_atop_mode, 
    dst_atop_mode, 
    xor_mode      
};

BlendProc findBlend(GBlendMode mode) {
    switch (mode) {
        case GBlendMode::kClear: 
            return clear_mode;

        case GBlendMode::kSrc:
            return src_mode;

        case GBlendMode::kDst:
            return dst_mode;

        case GBlendMode::kSrcOver:
            return src_over_mode;

        case GBlendMode::kDstOver:
            return dst_over_mode;

        case GBlendMode::kSrcIn:
            return src_in_mode;

        case GBlendMode::kDstIn:
            return dst_in_mode;

        case GBlendMode::kSrcOut:
            return src_out_mode;

        case GBlendMode::kDstOut:
            return dst_out_mode;

        case GBlendMode::kSrcATop:
            return src_atop_mode;

        case GBlendMode::kDstATop:
            return dst_atop_mode;

        case GBlendMode::kXor:
            return xor_mode;
    }
    return nullptr;
}

GPixel clear_mode(GPixel src, GPixel dst) {
    return GPixel_PackARGB(0, 0, 0, 0);
}

GPixel src_mode(GPixel src, GPixel dst) {
    return src;
}

GPixel dst_mode(GPixel src, GPixel dst) {
    return dst;
}

GPixel src_over_mode(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t dstA = GPixel_GetA(dst);

    if (srcA == 0) return dst;
    if (srcA == 255) return src;

    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    uint8_t invSrcA = 255 - srcA;
    uint8_t outA = std::min(255, srcA + ((dstA * invSrcA) >> 8));

    uint8_t outR = std::min(255, srcR + ((dstR * invSrcA) >> 8));
    uint8_t outG = std::min(255, srcG + ((dstG * invSrcA) >> 8));
    uint8_t outB = std::min(255, srcB + ((dstB * invSrcA) >> 8));

    return GPixel_PackARGB(outA, outR, outG, outB);
}

GPixel MyCanvas::blendPixel(GPixel dstPixel, GPixel srcPixel, GBlendMode blendMode) const {
    BlendProc blendProc = findBlend(blendMode);
    if (blendProc) {
        return blendProc(srcPixel, dstPixel);
    }
    return dstPixel;
}


GPixel dst_over_mode(GPixel src, GPixel dst) {
   uint8_t dstA = GPixel_GetA(dst);
    if (dstA == 0){
        return src;
    }
    if (dstA == 255){
        return dst;
    }

    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    uint8_t dstAlphaComp = 255 - dstA;

    uint8_t resultA = dstA + ((srcA * dstAlphaComp) >> 8);
    uint8_t resultR = dstR + ((srcR * dstAlphaComp) >> 8);
    uint8_t resultG = dstG + ((srcG * dstAlphaComp) >> 8);
    uint8_t resultB = dstB + ((srcB * dstAlphaComp) >> 8);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);

}

GPixel src_in_mode(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t dstA = GPixel_GetA(dst);

    if (srcA == 0 || dstA == 0) {
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);
    uint8_t resultA = (srcA * dstA) >> 8;
    uint8_t resultR = (srcR * dstA) >> 8;
    uint8_t resultG = (srcG * dstA) >> 8;
    uint8_t resultB = (srcB * dstA) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);

}

GPixel dst_in_mode(GPixel src, GPixel dst) {
     uint8_t srcA = GPixel_GetA(src);
    uint8_t dstA = GPixel_GetA(dst);

    if (srcA == 0 || dstA == 0) {
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);
    
    uint8_t resultA = (dstA * srcA) >> 8;
    uint8_t resultR = (dstR * srcA) >> 8;
    uint8_t resultG = (dstG * srcA) >> 8;
    uint8_t resultB = (dstB * srcA) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);

}

GPixel src_out_mode(GPixel src, GPixel dst) {
     uint8_t srcA = GPixel_GetA(src);
    if (srcA == 0){
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);
    uint8_t dstA = GPixel_GetA(dst);

    uint8_t dstAlphaComp = 255 - dstA;

    uint8_t resultA = (srcA * dstAlphaComp) >> 8;
    uint8_t resultR = (srcR * dstAlphaComp) >> 8;
    uint8_t resultG = (srcG * dstAlphaComp) >> 8;
    uint8_t resultB = (srcB * dstAlphaComp) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);

}

GPixel dst_out_mode(GPixel src, GPixel dst) {
    uint8_t dstA = GPixel_GetA(dst);
    if (dstA == 0){
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);
    uint8_t srcA = GPixel_GetA(src);

    uint8_t srcAlphaComp = 255 - srcA;

    uint8_t resultA = (dstA * srcAlphaComp) >> 8;
    uint8_t resultR = (dstR * srcAlphaComp) >> 8;
    uint8_t resultG = (dstG * srcAlphaComp) >> 8;
    uint8_t resultB = (dstB * srcAlphaComp) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

GPixel src_atop_mode(GPixel src, GPixel dst) {
       uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float dstAlpha = dstA / 255.0f;
    float srcAlphaComp = 1.0f - (srcA / 255.0f);

    uint8_t resultA = (int)(std::min(255.0f, (srcA * dstAlpha) + (dstA * srcAlphaComp) + 0.5f));
    uint8_t resultR = (int)((srcR * dstAlpha) + (dstR * srcAlphaComp) + 0.5f);
    uint8_t resultG = (int)((srcG * dstAlpha) + (dstG * srcAlphaComp) + 0.5f);
    uint8_t resultB = (int)((srcB * dstAlpha) + (dstB * srcAlphaComp) + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);

}

GPixel dst_atop_mode(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float srcAlpha = srcA / 255.0f;
    float dstAlphaComp = 1.0f - (dstA / 255.0f);

    uint8_t resultA = (int)(std::min(255.0f, (dstA * srcAlpha) + (srcA * dstAlphaComp) + 0.5f));
    uint8_t resultR = (int)((dstR * srcAlpha) + (srcR * dstAlphaComp) + 0.5f);
    uint8_t resultG = (int)((dstG * srcAlpha) + (srcG * dstAlphaComp) + 0.5f);
    uint8_t resultB = (int)((dstB * srcAlpha) + (srcB * dstAlphaComp) + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

GPixel xor_mode(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float srcAlphaComp = 1.0f - (srcA / 255.0f);
    float dstAlphaComp = 1.0f - (dstA / 255.0f);

    uint8_t resultA = (int)(srcA * dstAlphaComp + dstA * srcAlphaComp + 0.5f);
    uint8_t resultR = (int)(srcR * dstAlphaComp + dstR * srcAlphaComp + 0.5f);
    uint8_t resultG = (int)(srcG * dstAlphaComp + dstG * srcAlphaComp + 0.5f);
    uint8_t resultB = (int)(srcB * dstAlphaComp + dstB * srcAlphaComp + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

Edge clipEdge(const Edge& edge, int canvasWidth, int canvasHeight) {
    Edge clippedEdge = edge;

    if (clippedEdge.top < 0) {
        clippedEdge.xLeft = static_cast<int>(std::round(clippedEdge.m * -clippedEdge.top + clippedEdge.b));
        clippedEdge.top = 0;
    }

    if (clippedEdge.bottom > canvasHeight) {
        clippedEdge.bottom = canvasHeight;
    }

    if (clippedEdge.xLeft < 0) {
        clippedEdge.xLeft = 0;
        clippedEdge.top = static_cast<int>(std::round((0 - clippedEdge.b) / clippedEdge.m));
    }

    if (clippedEdge.xRight > canvasWidth) {
        clippedEdge.xRight = canvasWidth;
        clippedEdge.bottom = static_cast<int>(std::round((canvasWidth - clippedEdge.b) / clippedEdge.m));
    }

    return clippedEdge;
}


Edge makeEdge(const GPoint& p0, const GPoint& p1) {
    GPoint top = p0, bottom = p1;
    int winding = 1;
    if (p0.y > p1.y) {
        std::swap(top, bottom);
        winding = -1;
    }

    if (top.y == bottom.y) {
        return {}; // Degenerate edge
    }

    Edge edge;
    edge.m = (bottom.x - top.x) / (bottom.y - top.y);  // Slope
    edge.b = top.x - edge.m * top.y;  // X = mY + b, solve for b
    edge.top = static_cast<int>(std::round(top.y));
    edge.bottom = static_cast<int>(std::round(bottom.y));
    edge.currX = top.x;
    edge.winding = winding;

    return edge;
}
