// Eigen 3.4.0
// https://eigen.tuxfamily.org/index.php?title=Main_Page
#include "../External_libraries/Eigen/Dense"

using namespace Eigen;

struct Triangle
{
    Vector4f p[3];      // 3 points forming the triangle in model-space. Each point is (x, y, z, w).
};


Matrix4f createTranslationMatrix(Vector3f vector)
{
    Matrix4f matrix;
    matrix << 1, 0, 0, vector.x(),
              0, 1, 0, vector.y(),
              0, 0, 1, vector.z(),
              0, 0, 0,     1;

    return matrix;
}

// 3D-rotation, theta in radians
Matrix4f createRotationXMatrix(float theta)
{
    Matrix4f matrix;
    matrix << 1,     0,           0,        0,
              0, cosf(theta), -sinf(theta), 0,
              0, sinf(theta), cosf(theta),  0,
              0,     0,           0,        1;

    return matrix;     

}

// 3D-rotation, theta in radians
Matrix4f createRotationYMatrix(float theta)
{
    Matrix4f matrix;
    matrix << cosf(theta),  0, sinf(theta), 0,
                   0,       1,      0,      0,
              -sinf(theta), 0, cosf(theta), 0,
                   0,       0,      0,      1;

    return matrix;
}

// 3D-rotation, theta in radians
Matrix4f createRotationZMatrix(float theta)
{
    Matrix4f matrix;
    matrix << cosf(theta), -sinf(theta), 0, 0,
              sinf(theta), cosf(theta),  0, 0,
              0,           0,            1, 0,
              0,           0,            0, 1;

    return matrix;
}

Matrix4f createScaleMatrix(float scaleX, float scaleY, float scaleZ)
{
    Matrix4f matrix;
    matrix << scaleX, 0,      0,       0,
              0,      scaleY, 0,       0,
              0,      0,      scaleZ,  0,
              0,      0,      0,       1;

    return matrix;
}

/**
 * @brief Transformation to go from camera-space (3D) to screen-space (2D)
 * 
 * @param[in] fovDegree     Field of view in degrees
 * @param[in] screenWidth   Screen width in pixels
 * @param[in] screenHeight  Screen height in pixels
 * @param[in] nearPlaneZ    Nearest plane before which objects will be clipped
 * @param[in] farPlaneZ     not used for now...
 * @return Matrix4f 
 */
Matrix4f createProjectionMatrix(float fovDegree, float screenWidth, float screenHeight, float nearPlaneZ, float farPlaneZ)
{
    float aspectRatio = screenHeight / screenWidth;
    float fovRadian = fovDegree * 2 * M_PI / 360.0f;
    float screenDist = 1.0f / tanf(fovRadian / 2.0f);   // Distance at which the objects are projected on the screen
    Matrix4f matrix = Matrix4f::Zero();
    
    matrix(0,0) = aspectRatio * screenDist;
    matrix(1,1) = screenDist;
    matrix(2,2) = farPlaneZ / (farPlaneZ - nearPlaneZ);                     // Not used... make Z between 0 and farPlaneZ
    matrix(2,3) = (-farPlaneZ * nearPlaneZ) / (farPlaneZ - nearPlaneZ);     // Not used...
    matrix(3,2) = 1.0f;

    return matrix;
}


/**
 * @brief Tranformation to go from world-space to camera-space. The transformation
 *        is the inverse transformation of the placement of the camera in world-space.
 * 
 * @param[in] vEye    Where the camera is located in world-space
 * @param[in] vTarget Where we want to look at
 * @param[in] vUp     The up direction (+1 ou -1 vector)
 * @return The look at matrix
 */
Matrix4f createLookAtMatrix(const Vector3f& vEye, const Vector3f& vTarget, const Vector3f& vUp)
{
    // Calculate the new reference system (camera) in world-space
    Vector3f zAxis = (vTarget - vEye).normalized();       // Forward direction
    Vector3f xAxis = (vUp.cross(zAxis)).normalized();   // Right direction
    Vector3f yAxis = zAxis.cross(xAxis);           // up direction

    Matrix4f rotation
    {
       { xAxis.x(), yAxis.x(), zAxis.x(), 0.0f },
       { xAxis.y(), yAxis.y(), zAxis.y(), 0.0f },
       { xAxis.z(), yAxis.z(), zAxis.z(), 0.0f },
       {    0.0f,     0.0f,      0.0f,    1.0f }
    };

    // The inverse of the rotation matrix is the transpose because the matrix is orthonormalized
    Matrix4f rotationInverse = rotation.transpose();
    
    // The inverse of a translation is just the opposite translation
    // T(v)^-1 == T(-v)
    Matrix4f translationInverse = createTranslationMatrix({ -vEye.x(), -vEye.y(), -vEye.z() });
 
    // The normal operation is translation * rotation, but the inverse is
    // rotationInverse * translationInverse.
    return ( rotationInverse * translationInverse );

}

/**
 * @brief Find the intersection between a line and a plane that
 *        are known to intersect.
 * 
 * @param[in] planePoint    A point of the plane
 * @param[in] planeNormal   The normal of the plane (orientation)
 * @param[in] pointA        The inside point of the line
 * @param[in] pointB        The outside point of the line
 * @param[out] t            The fraction of (pointB - pointA) to go from pointA to the intersection
 * @return The point of intersection
 */
Vector3f lineIntersectPlane(const Vector3f& planePoint, const Vector3f& planeNormal, const Vector3f& pointA, const Vector3f& pointB, float& t)
{
    // Almost copy-pasted :)
    float plane_d = - planeNormal.dot(planePoint);      // No idea what plane_d represents
    float ad = pointA.dot(planeNormal);
    float bd = pointB.dot(planeNormal);
    t = (-plane_d - ad) / (bd - ad);        // The fraction of (pointB - pointA) to go from pointA to the intersection
    Vector3f vLine = pointB - pointA;
    Vector3f vPointAToIntersection = vLine * t;

    return pointA + vPointAToIntersection;
}


/**
 * @brief Return the triangles resulting from the clipping of a triangle against a plane
 * 
 * @param planePoint    A point of the plane
 * @param planeNormal   The normal of the plane
 * @param triangle      The triangle to clip
 * @return All the resulting triangles (between 0 and 2)
 */
std::vector<Triangle> clipTriangleAgainsPlane(Vector3f planePoint, Vector3f planeNormal, Triangle& triangle)
{
    std::vector<Triangle> clippedTriangles;
    planeNormal.normalized();   // Just to be sure

    // Shortest distance from point to plane including direction (+ -> inside, - -> outside)
    auto calulateDist = [&planePoint, &planeNormal](Vector3f point)
    {
        return (planeNormal.dot(point) - planeNormal.dot(planePoint));
    };

    // stupid, will see what is the performance impact...
    auto convert = [&triangle](Vector3f vertex3D)
    {
        Vector4f vertex4D;
        vertex4D.x() = vertex3D.x();
        vertex4D.y() = vertex3D.y();
        vertex4D.z() = vertex3D.z();
        vertex4D.w() = triangle.p[0].w();   // The transformations here don't change w
        return vertex4D;
    };

    // Use 3D coordinates instead
    Vector3f vertex[3];

    // Sort points if they are inside or outside the zone of interest
    // The points ordering is kept in each array.
    Vector3f* insidePoints[3]; uint32_t nbInsidePoints = 0u;
    Vector3f* outsidePoints[3]; uint32_t nbOutsidePoints = 0u;

    for (uint32_t i = 0u; i < 3; i++)
    {
        vertex[i] = triangle.p[i](seq(0,2));
        float dist = calulateDist(vertex[i]);
        if (dist >= 0.0f)
            insidePoints[nbInsidePoints++] = &vertex[i];
        else
            outsidePoints[nbOutsidePoints++] = &vertex[i];
    }


    // Case #1 : All points are outside the zone. The triangle is excluded.
    if (nbInsidePoints == 0)
    {}
    // Case #2 : All points are inside the zone. The triangle is kept.
    else if (nbInsidePoints == 3)
    {
        clippedTriangles.push_back(triangle);
    }
    // Case #3 : The triangle is clipped and result in 1 new triangle
    else if (nbInsidePoints == 1 && nbOutsidePoints == 2)
    {
        float t = 0.0f;
        Triangle newTriangle;

        //! Really important to keep the same points ordering as the initial triangle! //

        newTriangle.p[0] = convert(*insidePoints[0]);
        newTriangle.p[1] = convert(lineIntersectPlane(planePoint, planeNormal, *insidePoints[0], *outsidePoints[0], t));
        newTriangle.p[2] = convert(lineIntersectPlane(planePoint, planeNormal, *insidePoints[0], *outsidePoints[1], t));
        clippedTriangles.push_back(newTriangle);
    }
    // Case #4 : The triangle is clipped in 2 new triangles
    else if (nbInsidePoints == 2 && nbOutsidePoints == 1)
    {
        float t = 0.0f;
        Triangle newTriangle1, newTriangle2;

        //! Really important to keep the same points ordering as the initial triangle! //

        newTriangle1.p[0] = convert(*insidePoints[0]);
        newTriangle1.p[1] = convert(*insidePoints[1]);
        newTriangle1.p[2] = convert(lineIntersectPlane(planePoint, planeNormal, *insidePoints[0], *outsidePoints[0], t));

        newTriangle2.p[0] = convert(*insidePoints[1]);
        newTriangle2.p[1] = newTriangle1.p[2];      // Are p[1] and p[2] inverted ??
        newTriangle2.p[2] = convert(lineIntersectPlane(planePoint, planeNormal, *insidePoints[1], *outsidePoints[0], t));
        clippedTriangles.push_back(newTriangle1);
        clippedTriangles.push_back(newTriangle2);
    }

    return clippedTriangles;
}
