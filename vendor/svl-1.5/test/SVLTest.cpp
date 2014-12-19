/*
    File:           SVLTest.cpp

    Function:       Test program for the SVL lib.

    Author(s):      Andrew Willmott

    Copyright:      (c) 1995-2001, Andrew Willmott

*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "svl/SVL.h"

using namespace std;


// --- Prototypes -------------------------------------------------------------

Void TestBasics();

Void Test2DStuff();
Void Test3DStuff();
Void Test4DStuff();

Void TestHStuff2D();
Void TestHStuff3D();

Void TestND();
Void TestNDNumerical();

Void TestInit();
Void TestInput();

#define TYPE_NAME(X) type_name((X*)0)
inline char *type_name(Float*) { return("Float"); };
inline char *type_name(Double*) { return("Double"); };
inline char *type_name(Int*) { return("Int"); };


int main(int, char **)
{
    cout << "Testing SVL library, version " << SVL_VERSION << endl;
    cout << "Real type: " << TYPE_NAME(Real) << endl;

    cout << "----------------------------------" << endl;

//  TestInput();

    Test2DStuff();
    Test3DStuff();
    Test4DStuff();

    TestHStuff2D();
    TestHStuff3D();

    TestND();
    TestNDNumerical();

    TestInit();

    cout << "\n\n--- Finished! ---" << endl;

    return(0);
}


// --- Test routines ----------------------------------------------------------

Void TestHStuff2D()
{
    Vec2 x;
    Vec3 y;
    Mat2 M;

    cout << "\n+ TestHStuff2D\n\n";

    x = Vec2(1,2);
    cout << "x is: " << x << endl;
    M = Rot2(vl_halfPi);
    cout << "rot(pi/2) is: " <<  Rot2(vl_halfPi) << endl;
    x = Rot2(vl_halfPi) * x;
    cout << "x after rot(pi/2) is: " << x << endl;
    x = Scale2(Vec2(0.3, 0.2)) * x;
    cout << "x after scale(0.3, 0.2) is: " << x << endl;

    y = Vec3(x, 0.5);
    cout << "y is: " << y << endl;
    x = proj(y);
    cout << "proj(y) is: " << x << endl;

    x = proj(HRot3(1.3) * HTrans3(Vec2(1,1)) *
             HScale3(Vec2(1,2)) * Vec3(x, 1));
    cout << "HRot3(1.3) * HTrans3(Vec2(1,1)) * HScale3(Vec2(1,2)) * y = "
         << x << endl;
}

Void TestHStuff3D()
{
    Vec3 x;
    Vec4 y;

    cout << "\n+ TestHStuff3D\n\n";

    x = Vec3(1,2,3);

    cout << "rot(pi/2, vl_x) is: " <<  Rot3(vl_x, vl_halfPi) << endl;
    x = x * Rot3(vl_x, vl_halfPi);
    cout << "x after rot(pi/2, vl_x) is: " << x << endl;
    x = x * Scale3(Vec3(0.3, 0.2, 0.3));
    cout << "x after scale(0.3, 0.2, 0.3) is: " << x << endl;

    y = Vec4(x, 0.5);
    cout << "y is: " << y << endl;
    x = proj(y);
    cout << "proj(y) is: " << x << endl;

    x = proj(HRot4(vl_x, 1.3) * HTrans4(vl_1) * HScale4(Vec3(1,2,1)) * y);
    cout << "HRot4(vl_x, 1.3) * HTrans4(vl_1) "
        "* HScale4(Vec3(1,2,1)) * y = " << x;
}


Void Test2DStuff()
{
    Vec2 x(1,2);
    Vec2 y(5,6);

    cout << "\n+ Test2DStuff\n\n";

    cout << "x: " << x << ", y: " << y << "\n\n";

    cout << "x + y * (y * x * 2) : " << x + y * (y * x * 2) << endl;
    cout << "x dot y               : " << dot(x, y) << endl;

    cout << "cross(x)    : " << cross(x) << endl;
    cout << "len         : " << len(x) << endl;
    cout << "sqrlen      : " << sqrlen(x) << endl;
    cout << "norm        : " << norm(x) << endl;
    cout << "len of norm : " << len(norm(x)) << "\n\n";

    Mat2 M(1,2,3,4);
    Mat2 N; N.MakeDiag(2.0);

    cout << "M       : " << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << "\n" << endl;

    cout << "Vec2 consts: " << Vec2(vl_0) << Vec2(vl_x)
         << Vec2(vl_y) << Vec2(vl_1) << endl;
    cout << "Mat2 consts:\n" << Mat2(vl_Z) << endl << Mat2(vl_I)
         << endl << Mat2(vl_B) << "\n\n";

    M = Rot2(1.3) * Scale2(Vec2(2,1));

    cout << "M       : \n" << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << "\n" << endl;
}

Void Test3DStuff()
{
    Vec3 x(1,2,3);
    Vec3 y(5,6,7);

    cout << "\n+ Test3DStuff\n\n";

    cout << "x: " << x << ", y: " << y << "\n\n";

    cout << "x + y * (y * x * 2) : " << x + y * (y * x * 2) << endl;
    cout << "x dot y               : " << dot(x, y) << endl;

    cout << "cross(x,y)  : " << cross(x,y) << endl;
    cout << "cross(x, y) . x : " << dot(cross(x, y), x) << endl;
    cout << "cross(x, y) . y : " << dot(cross(x, y), y) << endl;
    cout << "len         : " << len(x) << endl;
    cout << "sqrlen      : " << sqrlen(x) << endl;
    cout << "norm        : " << norm(x) << endl;
    cout << "len of norm : " << len(norm(x)) << "\n\n";

    Mat3 M(1,2,3,3,2,1,2,1,3);
    Mat3 N; N.MakeDiag(2.0);

    cout << "M       : \n" << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << endl;

    cout << "Vec3 consts: " << Vec3(vl_0) << Vec3(vl_x)
         << Vec3(vl_y) << Vec3(vl_z) << Vec3(vl_1) << endl;
    cout << "Mat3 consts:\n" << Mat3(vl_Z) << endl << Mat3(vl_I)
         << endl << Mat3(vl_B) << "\n\n";

    M = Rot3(vl_y, 1.3) * Scale3(Vec3(2,4,2));

    cout << "M       :\n" << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << endl;
}

Void Test4DStuff()
{
    Vec4 x(1,2,3,4);
    Vec4 y(5,6,7,8);
    Vec4 z(1,0,0,0);

    cout << "\n+ Test4DStuff\n\n";

    cout << "x: " << x << ", y: " << y << ", z: " << z << "\n\n";

    cout << "x + y * (z * x * 2) : " << x + y * (z * x * 2) << endl;
    cout << "x dot y               : " << dot(x, y) << "\n\n";

    cout << "cross(x,y,z)   : " << cross(x,y,z) << endl;
    cout << "cross(x,y,z).x : " << dot(cross(x,y,z), x) << endl;
    cout << "cross(x,y,z).y : " << dot(cross(x,y,z), y) << endl;
    cout << "cross(x,y,z).z : " << dot(cross(x,y,z), z) << endl;
    cout << "len            : " << len(x) << endl;
    cout << "sqrlen         : " << sqrlen(x) << endl;
    cout << "norm           : " << norm(x) << endl;
    cout << "len of norm    : " << len(norm(x)) << "\n\n";


    Mat4 M(1,2,3,0, 2,3,0,5, 3,0,5,6, 0,5,6,7);
    Mat4 N; N.MakeBlock(2.0);

    cout << "M       : \n" << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << endl;

    cout << "Vec4 consts: " << Vec4(vl_0) << Vec4(vl_x) << Vec4(vl_y)
         << Vec4(vl_z) << Vec4(vl_w) << Vec4(vl_1) << endl;
    cout << "Mat4 consts:\n" << Mat4(vl_Z) << endl << Mat4(vl_I) << endl
         << Mat4(vl_B) << "\n\n";

    M = HScale4(Vec3(2,3,4));
    M *= HRot4(vl_y, 1.256);

    cout << "M       : " << M << endl;

    cout << "M * x   : " << M * x << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "adj     : " << adj(M) << endl;
    cout << "det     : " << det(M) << endl;
    cout << "trace   : " << trace(M) << endl;
    cout << "inv     : \n" << inv(M) << endl;
    cout << "M * inv : \n" << clamped(M * inv(M)) << endl;

}

Void TestND()
{
    cout << "\n+ TestND\n\n";

    Vec x(4, 1.0, 2.0, 3.0, 4.0);
    Vec sx(3, 1.0, 2.0, 3.0);
    Vec y(4, 5.0, 6.0, 7.0, 8.0);
    Vec z(4, 4.0, 3.0, 2.0, 1.0);

    Mat M(4, 3, 4.0, 3.0, 2.0, 1.0, 4.0, 3.0, 2.0, 1.0, 5.0,
           2.0, 1.0, 4.0, 3.0, 2.0, 1.0, 5.0);
    Mat N(4, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
           2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);

    cout << "M:\n" << M;
    cout << "\nN:\n" << N;
    cout << "\ntrans(N):\n" << trans(N);
    cout << "\nM + N:\n" << M + N;
    cout << "\nM * trans(N):\n" << M * trans(N);
    cout << "x: " << x << ", sx: " << sx << endl;

    z = x + y;
    cout << "z: " << z << endl;

    cout << "M * sx   : " << (M * sx) << endl;
    cout << "x * M   : " << x * M << endl;

    cout << "x: " << x << ", y: " << y << ", z: " << z << endl;

    cout << "x + y: " << x + y << endl;
    cout << "x * y: " << x * y << endl;
    cout << "x dot y: " << dot(x, y) << endl;
    cout << "x + (y dot z) * x * 2 : " << x + (dot(y, z) * x * 2.0) << endl;
    cout << "x + (y dot z) * x * 2 : " << ((float)2.0 * dot(y, z)) * x << endl;
    cout << "x + (y dot z) * x * 2 : " << x + dot(y, z) * x << endl;

    cout << "len : " << len(x) << endl;
    cout << "sqrlen : " << sqrlen(x) << endl;
    cout << "norm : " << norm(x) << endl;
    cout << "len of norm : " << len(norm(x)) << endl;
}

Void TestNDNumerical()
{
    Mat P(4, 4,
                1.0, 2.0, 3.0, 0.0,
                2.0, 3.0, 0.0, 5.0,
                3.0, 0.0, 5.0, 6.0,
                0.0, 5.0, 6.0, 7.0
            );
    Mat Q;
    Mat4    P1(
                1.0, 2.0, 3.0, 0.0,
                2.0, 3.0, 0.0, 5.0,
                3.0, 0.0, 5.0, 6.0,
                0.0, 5.0, 6.0, 7.0
            );
    Mat4    Q1;

    cout << "\n+ TestNDNumerical\n" << endl;

    cout << "P:\n";
    cout << P;

    cout << "\ninv(P):" << endl;
    Q = inv(P);
    cout << Q;

    cout << "\nP * inv(P):\n";
    cout << clamped(P * Q);

    cout << "\n\nReal:\n";

    cout << "P1: " << P1 << endl;
    cout << "\ninv(P1):\n";
    Q1 = inv(P1);
    cout << Q1 << endl;

    cout << "\nP1 * inv(P1): " << endl << clamped(P1 * Q1) << endl;
    cout << "\ninv(P) - inv(P1): " << endl << clamped(inv(P) - inv(P1));
    cout << endl << endl;
}

Void TestBasics()
{
    cout << "\n+ TestBasics\n\n";

    Assert(false, "Bogus assertion");
    Expect(false, "This is true");

    cout << "max is: " << Max(2.0, 5.0) << endl;
    cout << "min is: " << Min((Int) 2, (Int) 5) << endl;
}

Void TestInit()
{
    Vec v00(10, vl_zero);
    Vec v01(10, vl_one);
    Vec v02(10, vl_x);
    Vec v03(10, vl_axis(5));
    Vec v04(10); v04.MakeBlock(5.0);

    Mat m00(5, 5, vl_zero);
    Mat m01(5, 5, vl_one);
    Mat m02(5, 5); m02.MakeDiag(4.0);
    Mat m03(5, 5); m03.MakeBlock(2.2);

    cout << "--- Test init code -------------------------------------------"
         << endl << endl;

    cout << v00 << endl;
    cout << v01 << endl;
    cout << v02 << endl;
    cout << v03 << endl;
    cout << v04 << endl;

    cout << endl;

    cout << m00 << endl;
    cout << m01 << endl;
    cout << m02 << endl;
    cout << m03 << endl;
    cout << endl;

    cout << "testing comparisons..." << endl;
    cout << endl;
    cout << (Mat2(vl_0) == vl_0) << endl;
    cout << (Mat3(vl_0) == vl_0) << endl;
    cout << (Mat4(vl_0) == vl_0) << endl;
    cout << (Mat(8, 8, vl_0) == Mat(8, 8, vl_0)) << endl;

    cout << (Mat2(vl_1) == vl_0) << endl;
    cout << (Mat3(vl_1) == vl_0) << endl;
    cout << (Mat4(vl_1) == vl_0) << endl;
    cout << (Mat(8, 8, vl_1) == Mat(8, 8, vl_0)) << endl;
}

Void TestInput()
{
    Mat m;
    Vec v;

    while (1)
    {
        cout << "Input a vector: " << flush;
        if (cin >> v)
            cout << endl << v << endl << endl;
        else
        {
            cout << "failed." << endl;
            break;
        }

        cout << "Input a matrix: " << flush;
        if (cin >> m)
            cout << endl << m << endl;
        else
        {
            cout << "failed." << endl;
            break;
        }
    }

    exit (0);
}
