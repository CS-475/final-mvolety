#include "GFinal.h"
#include "include/GColor.h"
#include "include/GPoint.h"
#include "include/GPathBuilder.h"
#include "include/GMatrix.h"
#include "include/GShader.h"
#include "linear_utils.h"   
#include <iostream>
#include <memory>

class MyFinal : public GFinal {
public:
    std::shared_ptr<GShader> createVoronoiShader() override {
        GPoint p0 = {0,0}, p1 = {100,0};
        GColor colors[2] = { {1,1,0,0}, {1,0,1,0} }; 
        return GCreateLinearGradient(p0, p1, colors, 2, GTileMode::kClamp);
    }

    std::shared_ptr<GShader> createSweepGradient() override {
        GPoint p0 = {0,0}, p1 = {0,256};
        GColor colors[2] = { {1,0,1,0}, {1,1,1,1} }; 
        return GCreateLinearGradient(p0, p1, colors, 2, GTileMode::kClamp);
    }

    std::shared_ptr<GShader> createLinearPosGradient() override {
        GPoint p0 = {0,0}, p1 = {0,256};
        GColor colors[2] = { {1,0,1,0}, {1,1,1,1} }; 
        return GCreateLinearGradient(p0, p1, colors, 2, GTileMode::kClamp);
    }

    std::shared_ptr<GShader> createColorMatrixShader() override {
        GPoint p0 = {0, 0};
        GPoint p1 = {0, 256};
        GColor colors[2] = {
            {1, 0, 1, 0}, 
            {1, 1, 1, 1}  
        };
        return GCreateLinearGradient(p0, p1, colors, 2, GTileMode::kClamp);
    }

    void strokePolygon() override {
        GPathBuilder builder;
        builder.moveTo(10,10);
        builder.lineTo(100,10);
        builder.lineTo(100,100);
        builder.lineTo(10,100);
        builder.lineTo(10,10); 

        fStrokePath = builder.detach();
    }

    void drawQuadraticCoons() override {
        GPathBuilder builder;
        builder.moveTo(0,0);
        builder.lineTo(50,0);
        builder.lineTo(50,50);
        builder.lineTo(0,50);
        builder.lineTo(0,0);

        fCoonsPath = builder.detach(); 
    }

    private:
    std::shared_ptr<GPath> fStrokePath;
    std::shared_ptr<GPath> fCoonsPath;

};

std::unique_ptr<GFinal> GCreateFinal() {
    return std::make_unique<MyFinal>();
}
