#define OLC_PGE_APPLICATION
#define EIGEN_WARNINGS_DISABLED
#include "olcPixelGameEngine.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <deque>
#include <chrono>

using namespace Eigen;
using namespace olc;

/*
Author: Marc-Antoine Huet
Thanks to Javidx9 and others for olcPixelGameEngine.
Date : may 2022 

Notes:
- All Vector4f reprensent a point in 4D-space (x, y, z, w)
- All Vector3f reprensent a point in 3D-space (x, y, z) or a vector in 3D-space
- For simplicity, all matrix operations are done in 4D-space and all geometric operations are done in 3D-space.
*/


struct Triangle
{
    Vector4f p[3];      // 3 points forming the triangle in model-space. Each point is (x, y, z, w).
};


struct Mesh
{
    std::vector<Triangle> triangles;

    void createSimpleCubeModel()
    {
        Triangle triangle;

        // FRONT
        triangle.p[0] << 0.0f, 0.0f, 0.0f, 1.0f;  triangle.p[1] << 0.0f, 1.0f, 0.0f, 1.0f;  triangle.p[2] << 1.0f, 1.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 0.0f, 0.0f, 0.0f, 1.0f;  triangle.p[1] << 1.0f, 1.0f, 0.0f, 1.0f;  triangle.p[2] << 1.0f, 0.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        // RIGHT
        triangle.p[0] << 1.0f, 0.0f, 0.0f, 1.0f;  triangle.p[1] << 1.0f, 1.0f, 0.0f, 1.0f;  triangle.p[2] << 1.0f, 1.0f, 1.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 1.0f, 0.0f, 0.0f, 1.0f;  triangle.p[1] << 1.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 1.0f, 0.0f, 1.0f, 1.0f;
        triangles.push_back(triangle);

        // BACK
        triangle.p[0] << 1.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 1.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 0.0f, 1.0f, 1.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 1.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 0.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 0.0f, 0.0f, 1.0f, 1.0f;
        triangles.push_back(triangle);

        // LEFT
        triangle.p[0] << 0.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 0.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 0.0f, 1.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 0.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 0.0f, 1.0f, 0.0f, 1.0f;  triangle.p[2] << 0.0f, 0.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        // TOP
        triangle.p[0] << 0.0f, 1.0f, 0.0f, 1.0f;  triangle.p[1] << 0.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 1.0f, 1.0f, 1.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 0.0f, 1.0f, 0.0f, 1.0f;  triangle.p[1] << 1.0f, 1.0f, 1.0f, 1.0f;  triangle.p[2] << 1.0f, 1.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        // BOTTOM
        triangle.p[0] << 1.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 0.0f, 0.0f, 1.0f, 1.0f;  triangle.p[2] << 0.0f, 0.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);

        triangle.p[0] << 1.0f, 0.0f, 1.0f, 1.0f;  triangle.p[1] << 0.0f, 0.0f, 0.0f, 1.0f;  triangle.p[2] <<  1.0f, 0.0f, 0.0f, 1.0f;
        triangles.push_back(triangle);
    }

    /*  File requirements:
     *  Format : .obj
     *  Supported cathegories : v, f 
     *  The faces (f) can be 3-sided or 4-sides. If 4-sided, the polygon is split in 2 triangles.
    */
    bool loadFromObjectFile(std::string fileName, bool hasTexture = false)
    {
        std::ifstream file(fileName);
        if(!file.is_open())
        {
            std::cout << "Filepath error!\n";
            return false;
        }

        // Read the file and store the data temporary 
        std::vector<Vector4f> verticesList;
        std::vector<std::vector<int32_t>> facesList;

        std::string line;
        std::stringstream lineStream;
        uint32_t ite = 0u;
        while(getline(file, line))
        {
            lineStream.str("");
            lineStream.clear();

            if(line[0] == 'v' && line[1] == ' ')
            {
                lineStream << line;
                char junkChar;
                lineStream >> junkChar;

                Vector4f vertex;
                lineStream >> vertex.x() >> vertex.y() >> vertex.z();
                vertex.w() = 1.0f;
                verticesList.push_back(vertex);
            }
            else if(line[0] == 'f' && line[1] == ' ')
            {
                char junkChar;
                int32_t value;
                std::vector<int32_t> face;
                lineStream << line;
                lineStream >> junkChar;

                if(hasTexture)
                {
                    int32_t junkInt;
                    while(lineStream >> value)
                    {
                        face.push_back(value);
                        lineStream >> junkChar;
                        lineStream >> junkInt;
                    }
                }
                else
                {
                    while(lineStream >> value)
                    {
                        face.push_back(value);
                    }
                }

                // Split each 4-sided polygon in 2 triangles
                if(face.size() == 4)
                {
                    std::vector<int32_t> newFace { face[2], face[3], face[0] };
                    face.pop_back();    // Delete the last vertex of the face

                    facesList.push_back(face);
                    facesList.push_back(newFace);
                }
                else if(face.size() == 3)
                {
                    facesList.push_back(face);
                }
                else
                {
                    std::cout << "Number of vertices of a face is invalid\n";
                    return false;
                }
                
            }
            ite++;
        }

        // Create a triangle for each face with the corresponding vertex
        for (auto& face: facesList)
        {
            if (face.size() != 3)
            {
                std::cout << "The face doesn't have 3 vertices. Size = " << face.size() << std::endl;
                return false;
            }

            Triangle newTriangle { verticesList.at(face[0] - 1),     // index of face start at 1
                                   verticesList.at(face[1] - 1), 
                                   verticesList.at(face[2] - 1) };
            triangles.push_back(newTriangle); 
        }

        // Find z-min and z-max
        int32_t zMin = 0;
        int32_t zMax = 0;
        for(auto& vertex : verticesList)
        {
            if(vertex.z() > zMax)
                zMax = vertex.z();
            if(vertex.z() < zMin)
                zMin = vertex.z();
        }

        std::cout << "Number of polygons: " << triangles.size() << std::endl;
        std::cout << "zMin: " << zMin << "zMax: " << zMax << std::endl;

        return true;
    }
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
Vector3f lineIntersectPlane(Vector3f& planePoint, Vector3f& planeNormal, Vector3f& pointA, Vector3f& pointB, float& t)
{
    // Almost copy-pasted :)
    planeNormal.normalize();
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
Matrix4f createLookAtMatrix(Vector3f& vEye, Vector3f& vTarget, Vector3f& vUp)
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



class Engine3D : public olc::PixelGameEngine
{
public:
    Engine3D() {}

    void updateCamera(float elapsedTime);

    // Called at initialization
    bool OnUserCreate() override
    {
       
        // Cube: 12 polygons
        // SpaceShip: 106 polygons
        // Mountain: 5000 polygons
        // Hurricos: 11 800 polygons
        // Artisans Hub: 9800 polygons

        // model_.createSimpleCubeModel(); 
        std::string fileName;
        
        //fileName = "/Users/marc-antoinehuet/Desktop/Marc-Antoine/Programming projects/3D_engine/models/PlayStation - Spyro 2 Riptos Rage - Hurricos/Hurricos.obj";
        //fileName = "/Users/marc-antoinehuet/Desktop/Marc-Antoine/Programming projects/3D_engine/models/VideoShip.obj";
        //fileName = "/Users/marc-antoinehuet/Desktop/Marc-Antoine/Programming projects/3D_engine/models/Mountain.obj";
        fileName = "/Users/marc-antoinehuet/Desktop/Marc-Antoine/Programming projects/3D_engine/models/PlayStation - Spyro the Dragon - Artisans Hub/Artisans Hub.obj";

        if(!model_.loadFromObjectFile(fileName, true))
            return false;

        matProjection_ = createProjectionMatrix(90.0f, ScreenWidth(), ScreenHeight(), 0.5f, 2000.0f);
        
        // See justification where the transform is used
        formatForAPITransform_ = createScaleMatrix(static_cast<float>(ScreenWidth()), static_cast<float>(ScreenHeight()), 1.0f)  // Pixel coordinate
                                 * createScaleMatrix(0.5f, 0.5f, 1.0f)              // 0 to +2 ----> 0 to +1
                                 * createTranslationMatrix({ 1.0f, 1.0f, 0.0f })    // -1 to +1 ----> 0 to +2
                                 * createScaleMatrix(-1.0f, -1.0f, 1.0f);       // Invert coordinates

        float theta_ = 0.0f;
        float yaw_ = 0.0f;
        vCameraPos_ = { 0.0f, 0.0f, 0.0f };
        vLookDir_ = { 0.0f, 0.0f, 1.0f };

        return true;
    }

    // Called once per frame
    bool OnUserUpdate(float elapsedTime) override
    {
        auto start1 = std::chrono::steady_clock::now();
        updateCamera(elapsedTime);

        // Transformation from model-space to world-space
        Matrix4f modelToWorldTransform = createTranslationMatrix({ 0.0f, 0.0f, 2.0f})
                                         * createRotationXMatrix(-M_PI/2)
                                         * createRotationZMatrix(0) 
                                         * createRotationYMatrix(0)
                                         * createScaleMatrix(1.0f, 1.0f, 1.0f);

        std::vector<Triangle> trianglesToRaster;

        // Draw each triangle one at a time
        for (auto triangle : model_.triangles)
        {
            // Transformation from model-space to world-space
            for (uint32_t i = 0u; i < 3; i++)
            {
                triangle.p[i] = modelToWorldTransform * triangle.p[i];
            }

            // Only keep the triangle that are visible from the camera
            Vector3f vVertex1 = (triangle.p[1] - triangle.p[0])(seq(0,2));
            Vector3f vVertex2 = (triangle.p[2] - triangle.p[1])(seq(0,2));

            Vector3f vTriangleNormal = (vVertex1.cross(vVertex2)).normalized();
            Vector3f vOneTrianglePlanePoint = triangle.p[0](seq(0,2));
            Vector3f vCameraRay = vOneTrianglePlanePoint - vCameraPos_;     // Take any point of the triangle

            if (vCameraRay.dot(vTriangleNormal) < 0.0f)
            {
                // Could have been done with a tensor instead of a for loop, but...
                for (uint32_t i = 0u; i < 3; i++)
                {
                    // Tranformation from world-space to camera-space
                    triangle.p[i] = worldToCameraTransform_ * triangle.p[i];
                }

                std::vector<Triangle> clippedTriangles;
                clippedTriangles = clipTriangleAgainsPlane( {0.0f, 0.0f, 0.5f}, Vector3f::UnitZ(), triangle);
                for (auto& clippedTriangle : clippedTriangles)
                {
                    for (uint32_t i = 0u; i < 3; i++)
                    {
                        // Transformation from camera-space to screen-space (2D).
                        clippedTriangle.p[i] = (matProjection_ * clippedTriangle.p[i]);
                           
                        
                        // At this point, the coordinates are (x, y, z', 1) and the triangle is projected on 2D screen. X and Y are 
                        // between -1 and +1. But, the coordinates for the API are as follow:
                        //
                        // In world space:          API (drawing windoe):
                        //  Y                       |---- X
                        //  |                       |
                        //  |                       |
                        //  ----- X                 Y
                        //
                        // 1. The X and Y coordinates need to be inverted.
                        // 2. The X and Y coordinates need to be between +0 and +1.
                        // 3. The X and Y coordinates need to be in pixel values. X -> [0 - ScreenWidth]. Y -> [0 - ScreenHeight]

                        clippedTriangle.p[i] = formatForAPITransform_ * clippedTriangle.p[i];
                        clippedTriangle.p[i] /= clippedTriangle.p[i].w();
                    }

                    trianglesToRaster.push_back(clippedTriangle);
                }
            }
        }

        auto end1 = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for matrix operations: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count()
            << " us" << std::endl;

        auto start2 = std::chrono::steady_clock::now();
        // Clear the screen
        Clear(Pixel(0,0,0));

        // COULD BE DONE BEFORE formatForAPITransform !!

        // Sort the triangles based on the mean z to draw them in order.
        sort(trianglesToRaster.begin(), trianglesToRaster.end(), [](Triangle &t1, Triangle &t2)
		{
			float z1 = (t1.p[0].z() + t1.p[1].z() + t1.p[2].z()) / 3.0f;
			float z2 = (t2.p[0].z() + t2.p[1].z() + t2.p[2].z()) / 3.0f;
			return z1 > z2;
		});

        // Clip the triangles agains the for planes (X+, X-. Y+, Y-)
        for (auto& triangle : trianglesToRaster)
        {
            std::deque<Triangle> trianglesList;

            trianglesList.push_back(triangle);
            int32_t nbNewTriangles = 1;

            // Clip all the triangles in the queue against a plane and add all the new triangles to the queue.
            // When all the triangles of the queue have been processed, go to the new plane.
            for (int32_t planeId = 0; planeId < 4; planeId++)
            {
                while (nbNewTriangles > 0)
                {
                    Triangle test = trianglesList.front();
                    trianglesList.pop_front();
                    nbNewTriangles--;

                    std::vector<Triangle> clippedTriangles;

                    switch (planeId)
					{
                        case 0:	// Ymin
                            clippedTriangles = clipTriangleAgainsPlane({ 0.0f, 0.0f, 0.0f }, Vector3f::UnitY(), test); break;
                        case 1: // Ymax
                            clippedTriangles = clipTriangleAgainsPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, -1 * Vector3f::UnitY(), test); break;
                        case 2:	// Xmin
                            clippedTriangles = clipTriangleAgainsPlane({ 0.0f, 0.0f, 0.0f }, Vector3f::UnitX(), test); break;
                        case 3:	// Xmax
                            clippedTriangles = clipTriangleAgainsPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, -1 * Vector3f::UnitX(), test); break;
					}

                    for (auto& clippedTriangle : clippedTriangles)
                        trianglesList.push_back(clippedTriangle);
                }
                // All the triangles of the queue have been clipped agains the plane. Update the number of triangles to process for the new plane
                nbNewTriangles = trianglesList.size();
            }

            auto end2 = std::chrono::steady_clock::now();
            std::cout << "Elapsed time for geometric operations: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count()
            << " us" << std::endl;

            auto start3 = std::chrono::steady_clock::now();
            // Draw the triangles on the screen
            for (auto& triangleToRaster : trianglesList)
            {
                FillTriangle(triangleToRaster.p[0].x(), triangleToRaster.p[0].y(), 
                             triangleToRaster.p[1].x(), triangleToRaster.p[1].y(),  
                             triangleToRaster.p[2].x(), triangleToRaster.p[2].y(), Pixel(225, 25, 25));

                DrawTriangle(triangleToRaster.p[0].x(), triangleToRaster.p[0].y(), 
                             triangleToRaster.p[1].x(), triangleToRaster.p[1].y(),  
                             triangleToRaster.p[2].x(), triangleToRaster.p[2].y(), Pixel(255, 255, 255));

                std::stringstream textOverlay;
                textOverlay << "Rendering (ms) " << elapsedTime * 1000;
                DrawString(200, 0, textOverlay.str());
            }

            auto end3 = std::chrono::steady_clock::now();
            std::cout << "Elapsed time for draw: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count()
            << " us" << std::endl;
        }
        

        return true;
    }


private:
    Mesh model_;
    Matrix4f worldToCameraTransform_;       // Transformation from world-space to camera-space
    Matrix4f formatForAPITransform_;        // Use the same system of reference as the API
    Matrix4f matProjection_;                // Projection matrix
    float theta_;           // Transformation from model-space to world-space (in radians)
    float yaw_;             // Camera rotation in XZ plane (in radians)
    Vector3f vCameraPos_;   // Camera location in world-space
    Vector3f vLookDir_;     // Direction where we are looking at from the point of view of the camera
};


void Engine3D::updateCamera(float fElapsedTime)
{
    if (GetKey(Key::UP).bHeld)
        vCameraPos_.y() += 100.0f * fElapsedTime;	    // Move the camera up

    if (GetKey(Key::DOWN).bHeld)
        vCameraPos_.y() -= 100.0f * fElapsedTime;	    // Move the camera down

    // Maybe add left right movement ???

    Vector3f vForwardMove = vLookDir_ * 100.0f * fElapsedTime;  // Forward movement if a key was pressed

    if (GetKey(Key::W).bHeld)
        vCameraPos_ += vForwardMove;    // Move the camera forward

    if (GetKey(Key::S).bHeld)
        vCameraPos_ -= vForwardMove;    // Move the camera backward

    if (GetKey(Key::A).bHeld)
        yaw_ += 1.0f * fElapsedTime;    // Rotate the camera left

    if (GetKey(Key::D).bHeld)
        yaw_ -= 1.0f * fElapsedTime;    // Rotate the camera right

    // Transformation from world-space to camera-space
    Vector3f vUp = Vector3f::UnitY();          // The up direction never changes
    vLookDir_ = (createRotationYMatrix(yaw_) * Vector4f::UnitZ())(seq(0,2));    // Look in XZ plane
    Vector3f vTarget = vCameraPos_ + vLookDir_;

    worldToCameraTransform_ = createLookAtMatrix(vCameraPos_, vTarget, vUp);
}


// Used to mesure API performance
class RefGameEngine : public olc::PixelGameEngine
{
public:
    RefGameEngine()
    {}

public:
    bool OnUserCreate() override
    {
        // Called once at the start, so create things here
        std::string fileName = "/Users/marc-antoinehuet/Desktop/Marc-Antoine/Programming projects/3D_engine/models/PlayStation - Spyro the Dragon - Artisans Hub/High.png";
        sprite_ = new Sprite();
        if(sprite_->LoadFromFile(fileName) != 1)
        {
            std::cout << "Sprite file load error!\n";
            return false;
        }

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {
        auto start = std::chrono::steady_clock::now();

        // This routine takes ~10-11 ms for 800x800
        for (int i = 0; i < 800; i++)
        {
            for (int j = 0; j < 800; j++)
            {
                Draw(i, j, sprite_->GetPixel(250, 250));
            }
        }

        auto end = std::chrono::steady_clock::now();

        std::cout << "Elapsed time in microseconds: "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << " us" << std::endl;

        return true;
    }

private:
    Sprite* sprite_;
};


int main(int argc, char const* argv[])
{
    Engine3D engine3D;
    if (engine3D.Construct(800, 800, 1, 1))
        engine3D.Start();

    // RefGameEngine refGameEngine;
    // if (refGameEngine.Construct(800, 800, 1, 1))
    //     refGameEngine.Start();

    return 0;
}
