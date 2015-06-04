#pragma once

#include "vector3.h"

typedef struct Matrix3 {
    double x[3][3];
} Matrix3;

void m3_new(Matrix3 *const m, double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
    m->x[0][0] = xx;
    m->x[1][0] = xy;
    m->x[2][0] = xz;

    m->x[0][1] = yx;
    m->x[1][1] = yy;
    m->x[2][1] = yz;

    m->x[0][2] = zx;
    m->x[1][2] = zy;
    m->x[2][2] = zz;
}

Matrix3 m3_transpose(const Matrix3 *const m) {
    Matrix3 mt;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            mt.x[i][j] = m->x[j][i];
        }
    }
    return mt;
}

double m3_det(const Matrix3 *const m) {
    double a0 = m->x[1][1]*m->x[2][2] - m->x[1][2]*m->x[2][1];
    double a1 = m->x[1][2]*m->x[2][0] - m->x[1][0]*m->x[2][2];
    double a2 = m->x[1][0]*m->x[2][1] - m->x[1][1]*m->x[2][0];
    return m->x[0][0]*a0 + m->x[0][1]*a1 + m->x[0][2]*a2;
}

Matrix3 m3_inv(const Matrix3 *const m) {
    double det = 1.0/m3_det(m);
    Matrix3 inv;
    inv.x[0][0] = det*(m->x[1][1]*m->x[2][2] - m->x[1][2]*m->x[2][1]);
    inv.x[0][1] = det*(m->x[1][2]*m->x[2][0] - m->x[1][0]*m->x[2][2]);
    inv.x[0][2] = det*(m->x[1][0]*m->x[2][1] - m->x[1][1]*m->x[2][0]);

    inv.x[1][0] = det*(m->x[0][1]*m->x[2][2] - m->x[0][2]*m->x[2][1]);
    inv.x[1][1] = det*(m->x[0][2]*m->x[2][0] - m->x[0][0]*m->x[2][2]);
    inv.x[1][2] = det*(m->x[0][0]*m->x[2][1] - m->x[0][1]*m->x[2][0]);

    inv.x[2][0] = det*(m->x[0][1]*m->x[1][2] - m->x[0][2]*m->x[1][1]);
    inv.x[2][1] = det*(m->x[0][2]*m->x[1][0] - m->x[0][0]*m->x[1][2]);
    inv.x[2][2] = det*(m->x[0][0]*m->x[1][1] - m->x[0][1]*m->x[1][0]);

    return inv;
}

Vector3 m3_mv(const Matrix3 *const m, const Vector3 *const v) {
    Vector3 mv;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            mv.x[i] = m->x[j][i] * v->x[j];
        }
    }
    return mv;
}

Matrix3 m3_mm(const Matrix3 *const m1, const Matrix3 *const m2) {
    Matrix3 mm;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            mm.x[i][j] = 0.0;
            for (int k=0; k<3; k++) {
                mm.x[i][j] = m1->x[i][k] * m2->x[k][j];
            }
        }
    }
    return mm;
}

bool m3_eq(const Matrix3 *const m1, const Matrix3 *const m2, double tol) {
    bool eq = true;
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            eq = abs(m1->x[i][j] - m1->x[i][j]) < tol;
            if (!eq) break;
        }
    }
    return eq;
}
