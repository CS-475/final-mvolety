#ifndef UTILS_H
#define UTILS_H

#include "include/GPixel.h"
#include "include/GBlendMode.h"
#include "include/GPoint.h" 
#include <cmath>
#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include <stack>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) , fRect(GIRect::WH(device.width(), device.height())){
        fMatrixStack.push(GMatrix());
    }

    void clear(const GColor& color) override;\
    void drawRect(const GRect& rect, const GPaint& paint) override;
    void drawConvexPolygon(const GPoint[], int count, const GPaint&) override;
    void drawPath(const GPath& path, const GPaint& paint) override;
    void save() override;
    GPixel blendPixel(GPixel dstPixel, GPixel srcPixel, GBlendMode blendMode) const;
    void restore() override;
    void concat(const GMatrix&) override;
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
                  int count, const int indices[], const GPaint& paint) override;
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                  int level, const GPaint& paint) override;

private:
  void drawTriangle(const GPoint& p0, const GPoint& p1, const GPoint& p2,
                      const GColor* c0, const GColor* c1, const GColor* c2,
                      const GPoint* t0, const GPoint* t1, const GPoint* t2,
                      const GPaint& paint);
     void computeBarycentric(int x, int y, const GPoint& p0, const GPoint& p1, const GPoint& p2, float& alpha, float& beta, float& gamma);
    const GBitmap fDevice;
    GMatrix fMatrix;
    std::stack<GMatrix> fMatrixStack;
    const GIRect fRect;
};

struct Edge {
    float m;    
    float b;    
    int top;    
    int bottom; 
    float xLeft;
    int xRight; 
    int currX; 
    int winding;

    bool isValid(int y) const {
        return (y >= top && y < bottom);  
    }

    float computeX(int y) const {
        return m * y + b;  
    }

    bool isUseful() const {
        return top < bottom;
    }
};

Edge clipEdge(const Edge& edge, int canvasWidth, int canvasHeight);
Edge makeEdge(const GPoint& p0, const GPoint& p1);
GPixel convertColorToPix(const GColor& color);
GPixel clear_mode(GPixel src, GPixel dst);
GPixel src_mode(GPixel src, GPixel dst);
GPixel dst_mode(GPixel src, GPixel dst);
GPixel src_over_mode(GPixel src, GPixel dst);
GPixel dst_over_mode(GPixel src, GPixel dst);
GPixel src_in_mode(GPixel src, GPixel dst);
GPixel dst_in_mode(GPixel src, GPixel dst);
GPixel src_out_mode(GPixel src, GPixel dst);
GPixel dst_out_mode(GPixel src, GPixel dst);
GPixel src_atop_mode(GPixel src, GPixel dst);
GPixel dst_atop_mode(GPixel src, GPixel dst);
GPixel xor_mode(GPixel src, GPixel dst);

using BlendProc = GPixel (*)(GPixel src, GPixel dst);
BlendProc findBlend(GBlendMode mode);

#endif