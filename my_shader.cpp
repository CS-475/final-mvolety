#include "shader_utils.h"
#include <algorithm> 
#include <cmath>      

MyShader::MyShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tileMode)
    : fDevice(bitmap), fLocalMatrix(localMatrix), fTileMode(tileMode) {}

bool MyShader::isOpaque() {
    for (int y = 0; y < fDevice.height(); ++y) {
        for (int x = 0; x < fDevice.width(); ++x) {
            if (GPixel_GetA(*fDevice.getAddr(x, y)) != 255) {
                return false;
            }
        }
    }
    return true;
}


bool MyShader::setContext(const GMatrix& matrix) {
    GMatrix combined = GMatrix::Concat(matrix, fLocalMatrix);
    fInverse = combined.invert();
    return fInverse.has_value(); 
}


void MyShader::shadeRow(int x, int y, int count, GPixel row[]) {
    GPoint srcPoint;
    GPoint dstPoint;

    for (int i = 0; i < count; ++i) {
        dstPoint = {x + i + 0.5f, y + 0.5f};
        fInverse->mapPoints(&srcPoint, &dstPoint, 1);

        int srcX = applyTileMode(srcPoint.x, fDevice.width());
        int srcY = applyTileMode(srcPoint.y, fDevice.height());

        if (srcX < 0 || srcX >= fDevice.width() || srcY < 0 || srcY >= fDevice.height()) {
            row[i] = GPixel_PackARGB(0, 0, 0, 0); 
        } else {
            row[i] = *fDevice.getAddr(srcX, srcY);
        }
    }
}


int MyShader::applyTileMode(float coord, int size) const {
    if (size <= 0) return 0; 

    float t = coord / size;
    switch (fTileMode) {
        case GTileMode::kClamp:
            t = std::max(0.0f, std::min(1.0f, t));
            break;
        case GTileMode::kRepeat:
            t = t - std::floor(t);  
            break;
        case GTileMode::kMirror:
            t = mirrored_t(t); 
            break;
    }
    return static_cast<int>(t * size) % size; 
}


float MyShader::mirrored_t(float t) const {
    int n = static_cast<int>(std::floor(t));
    t = t - n;
    if (n % 2 == 1) {
        t = 1.0f - t;
    }
    return t;
}

std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode mode) {
    if (bitmap.width() <= 0 || bitmap.height() <= 0) {
        return nullptr;
    }
    return std::make_shared<MyShader>(bitmap, localMatrix, mode);
}
