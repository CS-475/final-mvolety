#include "include/GPathBuilder.h"
#include "include/GMatrix.h"

void GPathBuilder::addRect(const GRect& rect, GPathDirection dir) {
    this->moveTo(rect.left, rect.top);  
    if (dir == GPathDirection::kCW) {
        this->lineTo(rect.right, rect.top);    
        this->lineTo(rect.right, rect.bottom); 
        this->lineTo(rect.left, rect.bottom);  
    } else { 
        this->lineTo(rect.left, rect.bottom);  
        this->lineTo(rect.right, rect.bottom); 
        this->lineTo(rect.right, rect.top);    
    }
    this->lineTo(rect.left, rect.top);  
}

void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    if (count < 2) {
        return; 
    }
    this->moveTo(pts[0]);  
    for (int i = 1; i < count; ++i) {
        this->lineTo(pts[i]);  
    }
}

inline void GPathBuilder::transform(const GMatrix& matrix) {
    for (size_t i = 0; i < fPts.size(); ++i) {
        fPts[i] = matrix * fPts[i];
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection dir) {
    const float kSqrt2Over2 = static_cast<float>(std::sqrt(2.0));
    float r = radius;
    
    GPoint points[8] = {
        {center.x + r, center.y},                      
        {center.x + r * kSqrt2Over2, center.y + r * kSqrt2Over2}, 
        {center.x, center.y + r},                     
        {center.x - r * kSqrt2Over2, center.y + r * kSqrt2Over2},  
        {center.x - r, center.y},                        
        {center.x - r * kSqrt2Over2, center.y - r * kSqrt2Over2},  
        {center.x, center.y - r},                     
        {center.x + r * kSqrt2Over2, center.y - r * kSqrt2Over2},  
    };

    moveTo(points[0]);
    for (int i = 1; i < 8; i += 2) {
        lineTo(points[i]);
        lineTo(points[(i + 1) % 8]);  
    }
    lineTo(points[0]);
}