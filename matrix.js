//
// Vector and Matrix routines for client side WebGL Math.
// Author: Wayne O. Cochran wcochran@vancouver.wsu.edu
//

function dot3(U, V) {
    return U[0]*V[0] + U[1]*V[1] + U[2]*V[2];
}

function cross3(U, V) {
    var UxV = [];
    UxV[0] = U[1]*V[2] - U[2]*V[1];
    UxV[1] = U[2]*V[0] - U[0]*V[2];
    UxV[2] = U[0]*V[1] - U[1]*V[0];
    return UxV;
}

function scale3(s, V) {
    V[0] *= s; V[1] *= s; V[2] *= s; 
    return V;
}

function norm3(V) {
    var s = 1.0/Math.sqrt(dot3(V,V));
    return scale3(s,V);
}

//
// Represents 4x4 WebGL Transformation on client.
// GL expects matrices to be stored in column-major order.
// Use the 'array' property to set a GLSL mat4 uniform var; e.g.:
//    gl.uniformMatrix4fv(program.ModelViewProj, false, MVP.array);
//
function Matrix4x4() {
    this.array = [1, 0, 0, 0,   // stored in column-major order
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1];
}

Matrix4x4.prototype.elem = function(row, col) {  // column-major order index
    return this.array[row + 4*col];
}

Matrix4x4.prototype.setElem = function(row, col, val) {
    this.array[row + 4*col] = val;
    return val;
}

Matrix4x4.prototype.copy = function(B) {
    for (var i = 0; i < 4*4; i++)
        this.array[i] = B.array[i];
    return this;
}

Matrix4x4.prototype.identity = function() {
    for (var r = 0; r < 4; r++)
        for (var c = 0; c < 4; c++)
            this.setElem(r,c, (r == c) ? 1 : 0);
    return this;
}

Matrix4x4.prototype.set4x4 = function(a00, a01, a02, a03,
                                      a10, a11, a12, a13,
                                      a20, a21, a22, a23,
                                      a30, a31, a32, a33) {
    this.setElem(0,0,a00); 
    this.setElem(0,1,a01); 
    this.setElem(0,2,a02); 
    this.setElem(0,3,a03);
    this.setElem(1,0,a10);
    this.setElem(1,1,a11);
    this.setElem(1,2,a12);
    this.setElem(1,3,a13);
    this.setElem(2,0,a20);
    this.setElem(2,1,a21);
    this.setElem(2,2,a22);
    this.setElem(2,3,a23);
    this.setElem(3,0,a30);
    this.setElem(3,1,a31);
    this.setElem(3,2,a32);
    this.setElem(3,3,a33);
    return this;
}

Matrix4x4.prototype.set3x4 = function(a00, a01, a02, a03,
                                      a10, a11, a12, a13,
                                      a20, a21, a22, a23) {
    this.set4x4(a00, a01, a02, a03,
                a10, a11, a12, a13,
                a20, a21, a22, a23,
                0,   0,   0,   1);
    return this;
}               

Matrix4x4.prototype.set3x3 = function(a00, a01, a02,
                                      a10, a11, a12,
                                      a20, a21, a22) {
    this.set4x4(a00, a01, a02, 0,
                a10, a11, a12, 0,
                a20, a21, a22, 0,
                0,   0,   0,   1);
    return this;
}               

Matrix4x4.prototype.mult = function(B) {
    var AB = new Matrix4x4;
    for (var r = 0; r < 4; r++)
        for (var c = 0; c < 4; c++) {
            var s = 0;
            for (var i = 0; i < 4; i++)
                s += this.elem(r,i)*B.elem(i,c);
            AB.setElem(r,c,s);
        }
    return AB;
}

//
//  A <-- A*B
//
Matrix4x4.prototype.concat = function(B) {
    var AB = this.mult(B);
    return this.copy(AB);
}

//
// Scale Transformation.
// M <-- M * S(sx,sy,sz)
//
Matrix4x4.prototype.scale = function(sx, sy, sz) {
    S = new Matrix4x4;
    S.set3x3(sx, 0,  0,
             0,  sy, 0,
             0,  0,  sz);
    return this.concat(S);
}

//
// Translation Transformation.
// M <-- M * T(dx,dy,dz)
//
Matrix4x4.prototype.translate = function(dx, dy, dz) {
    T = new Matrix4x4;
    T.set3x4(1, 0, 0, dx,
             0, 1, 0, dy, 
             0, 0, 1, dz);
    return this.concat(T);
}

//
// Rotate about vector(x,y,z) Transformation.
// M <-- M * R(angle, x,y,z)
//
Matrix4x4.prototype.rotate = function(angle_degrees, x, y, z) {
    var p = 1/Math.sqrt(x*x + y*y + z*z);
    x *= p; y *= p; z *= p;
    var angle = angle_degrees * (Math.PI/180);
    var c = Math.cos(angle);
    var s = Math.sin(angle);
    var c_ = 1 - c;
    var zc_ = z*c_;
    var yc_ = y*c_;
    var xzc_ = x*zc_;
    var xyc_ = x*y*c_;
    var yzc_ = y*zc_;
    var xs = x*s;
    var ys = y*s;
    var zs = z*s;
    var S = new Matrix4x4;
    S.set3x3(x*x*c_ + c,  xyc_ - zs,   xzc_ + ys,
             xyc_ + zs,   y*yc_ + c,   yzc_ - xs,
             xzc_ - ys,   yzc_ + xs,   z*zc_ + c);
    return this.concat(S);
}

//
// "Look At" Transformation.
// Typically used to initialiaze View Transformation Matrix.
// M <-- M * Lookat
//   eye : position of camera
//   center : point camera is looking at
//   up : general "up" direction
//
Matrix4x4.prototype.lookat = function(eyex, eyey, eyez,
                                      centerx, centery, centerz,
                                      upx, upy, upz) {
    var F = [centerx - eyex, centery - eyey, centerz - eyez];
    norm3(F);
    var U = [upx, upy, upz];
    var S = cross3(F,U);
    norm3(S);
    U = cross3(S,F);
    var R = new Matrix4x4;
    R.set3x3( S[0],  S[1],  S[2],
              U[0],  U[1],  U[2],
             -F[0], -F[1], -F[2]);
    this.concat(R);
    return this.translate(-eyex, -eyey, -eyez);
}

//
// Orthonormal Projection
// M <-- M * Ortho
//
Matrix4x4.prototype.ortho = function(left, right,
                                     bottom, top,
                                     hither, yon) {
    var S = new Matrix4x4;
    S.set3x4(2/(right - left), 0,  0, -(right + left)/(right - left),
             0, 2/(top - bottom), 0, -(top + bottom)/(top - bottom),
             0, 0, -2/(yon-hither), -(yon + hither)/(yon - hither));
    return this.concat(S);
}

//
// 2D Projection Transformation
// M <-- M * Ortho2D
// visible z-values in [-1,+1].
//
Matrix4x4.prototype.ortho2D = function(left, right,
                                       bottom, top) {
    return this.ortho(M, left, right, bottom, top, -1, +1);
}

//
// Simple Perperspective Projection Trnaformation
// M <-- M * Persp
//  fovy : field of view angle along y-axis
//  aspect : aspect ratio of projection plane
//  zNear, zFar : near and far clipping planes
//
Matrix4x4.prototype.perspective = function(fovy_degrees, aspect, zNear, zFar) {
    var fovy = fovy_degrees*(Math.PI/180);
    var f = 1/Math.tan(fovy/2);
    var s = 1/(zNear - zFar);
    P = new Matrix4x4;
    P.set4x4(f/aspect, 0, 0, 0,
             0, f, 0, 0,
             0, 0, (zFar + zNear)*s, 2*zFar*zNear*s,
             0, 0, -1, 0);
    return this.concat(P);
}

//
// Normal Transformation
// Used to transform surface normals.
// Computes inverse-transpose of upper 3x3 of M.
// Returns array of 9 values representing 3x3 matrix 
// stored in column-major order.
//   XXX place GLSL set uniform example here
//
Matrix4x4.prototype.normal = function() {
    var M = this;
    var determinant =    
        +M.elem(0,0)*(M.elem(1,1)*M.elem(2,2) - M.elem(2,1)*M.elem(1,2))
        -M.elem(0,1)*(M.elem(1,0)*M.elem(2,2) - M.elem(1,2)*M.elem(2,0))
        +M.elem(0,2)*(M.elem(1,0)*M.elem(2,1) - M.elem(1,1)*M.elem(2,0));
    var invDet = 1.0/determinant;
    var normalMatrix = [];
    var N = function(row,col,val) { normalMatrix[col*3 + row] = val; }
    N(0,0, (M.elem(1,1)*M.elem(2,2) - M.elem(2,1)*M.elem(1,2))*invDet);
    N(1,0,-(M.elem(0,1)*M.elem(2,2) - M.elem(0,2)*M.elem(2,1))*invDet);
    N(2,0, (M.elem(0,1)*M.elem(1,2) - M.elem(0,2)*M.elem(1,1))*invDet);
    N(0,1,-(M.elem(1,0)*M.elem(2,2) - M.elem(1,2)*M.elem(2,0))*invDet);
    N(1,1, (M.elem(0,0)*M.elem(2,2) - M.elem(0,2)*M.elem(2,0))*invDet);
    N(2,1,-(M.elem(0,0)*M.elem(1,2) - M.elem(1,0)*M.elem(0,2))*invDet);
    N(0,2, (M.elem(1,0)*M.elem(2,1) - M.elem(2,0)*M.elem(1,1))*invDet);
    N(1,2,-(M.elem(0,0)*M.elem(2,1) - M.elem(2,0)*M.elem(0,1))*invDet);
    N(2,2, (M.elem(0,0)*M.elem(1,1) - M.elem(1,0)*M.elem(0,1))*invDet);
    return normalMatrix;
}

//
// Reflection Transformation.
// M <-- M * Reflect
//   reflection plane Ax + By + Cd + Z = 0 
//
Matrix4x4.prototype.reflect = function(plane) {
    var R = new Matrix3x3;
    var A = plane[0];
    var B = plane[1];
    var C = plane[2];
    var D = plane[3];
    R.set3x4(1-2*A*A,  -2*A*B,  -2*A*C, -2*A*D,
              -2*B*A, 1-2*B*B,  -2*B*C, -2*B*D,
              -2*C*A,  -2*C*B, 1-2*C*C, -2*C*D);
    return this.concat(R);
}

//
// Shadow transformation
// M <-- M * Shadow
//  light : homogeneous position of light source
//  plane : plane shadow is projected onto.
//
Matrix4x4.prototype.shadow = function(light, plane) {
    var dot =
        plane[0]*light[0] + plane[1]*light[1] +
        plane[2]*light[2] + plane[3]*light[3];
    
    var S = new Matrix4x4;
    
    S.setElem(0,0,dot - light[0]*plane[0]);
    S.setElem(0,1,    - light[0]*plane[1]);
    S.setElem(0,2,    - light[0]*plane[2]);
    S.setElem(0,3,    - light[0]*plane[3]);
        
    S.setElem(1,0,    - light[1]*plane[0]);
    S.setElem(1,1,dot - light[1]*plane[1]);
    S.setElem(1,2,    - light[1]*plane[2]);
    S.setElem(1,3,    - light[1]*plane[3]);
    
    S.setElem(2,0,    - light[2]*plane[0]);
    S.setElem(2,1,    - light[2]*plane[1]);
    S.setElem(2,2,dot - light[2]*plane[2]);
    S.setElem(2,3,    - light[2]*plane[3]);
        
    S.setElem(3,0,    - light[3]*plane[0]);
    S.setElem(3,1,    - light[3]*plane[1]);
    S.setElem(3,2,    - light[3]*plane[2]);
    S.setElem(3,3,dot - light[3]*plane[3]);
    
    return this.concat(S);
}

//
// Stack class for saving / restoring transformation matrics.
//
var Matrix4x4Stack = function() {
    this.stack = [];
}

Matrix4x4Stack.prototype.push = function(M) {
    var C = new Matrix4x4;
    C.copy(M);
    this.stack.push(C);
}

Matrix4x4Stack.prototype.pop = function(M) {
    var C = this.stack.pop();  
    M.copy(C);
}

