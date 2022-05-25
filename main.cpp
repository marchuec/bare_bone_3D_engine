#define OLC_PGE_APPLICATION
#define EIGEN_WARNINGS_DISABLED
#include "olcPixelGameEngine.h"
#include <iostream>
#include <fstream>
#include <deque>
#include <chrono>

// Matrix library and custom matrices
#include "MatricesAndUtilities.h"

using namespace olc;

/*
Author: Marc-Antoine Huet
Thanks to Javidx9 and others for olcPixelGameEngine.
Date : may 2022 

Notes:
- All Vector4f reprensent a point in 4D-space (x, y, z, w)
- All Vector3f reprensent a point in 3D-space (x, y, z) or a vector in 3D-space
- For simplicity, all matrix operations are done in 4D-space and all geometric operations are done in 3D-space.

- To fix: need to change the line  < #include "../External_libraries/libpng16/png.h" > (line 4084) in oldPixelGameEngine.h to your own libpng16 installation !!
*/


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

        std::cout << "** Model usefull data ** \n"
                  << "Number of vertices: " << verticesList.size() << std::endl
                  << "Number of polygons: " << triangles.size() << std::endl
                  << "zMin: " << zMin << "zMax: " << zMax << std::endl;

        return true;
    }
};


class Engine3D : public olc::PixelGameEngine
{
public:
    Engine3D() {}

    // Called at initialization
    bool OnUserCreate() override
    {
        // All the models available
        enum class Model { artisans_hub, hurricos, cube, videoShip, mountain, axis };

        /*
        Cube: 36 vertices / 12 polygons
        SpaceShip: 55 vertices / 106 polygons
        Mountain: 2441 vertices / 5000 polygons
        Hurricos: 5150 vertices / 11 800 polygons
        Artisans Hub: 6230 vertices / 9800 polygons 
        */

        Model modelToLoad = Model::artisans_hub;

        bool status = false;
        switch (modelToLoad)
        {
            case Model::artisans_hub :
                status = model_.loadFromObjectFile("../models/PlayStation - Spyro the Dragon - Artisans Hub/Artisans Hub.obj", true);
                break;
            case Model::hurricos :
                status = model_.loadFromObjectFile("../models/PlayStation - Spyro 2 Riptos Rage - Hurricos/Hurricos.obj", true);
                break;
            case Model::cube : 
                model_.createSimpleCubeModel(); status = true;
                break;
            case Model::mountain :
                status = model_.loadFromObjectFile("../models/Mountain.obj", false);
                break;
            case Model::axis :
                status = model_.loadFromObjectFile("../models/Axis.obj", false);
                break;
            case Model::videoShip :
                status = model_.loadFromObjectFile("../models/VideoShip.obj", false);
                break;
            default:
                return false;
        }

        if(!status)
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
        elapsedTime_ = elapsedTime;

        // Apply the matrix operations to all the triangles of the model
        std::vector<Triangle> result = applyMatrixOps();

        // Clip all the triangles agains the 
        clipAndDraw(result);

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
    float elapsedTime_;     // Elapsed time since last frame

    void updateCamera(void);
    void clipAndDraw(std::vector<Triangle>& triangles);
    std::vector<Triangle> applyMatrixOps(void);
};


std::vector<Triangle> Engine3D::applyMatrixOps(void)
{
    auto start1 = std::chrono::steady_clock::now();
    updateCamera();

    // Transformation from model-space to world-space
    Matrix4f modelToWorldTransform = createTranslationMatrix({ 0.0f, 0.0f, 300.0f})
                                        * createRotationXMatrix(-M_PI/2)
                                        * createRotationZMatrix(0) 
                                        * createRotationYMatrix(0)
                                        * createScaleMatrix(1.0f, 1.0f, 1.0f);

    std::vector<Triangle> processedTriangles;

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

                processedTriangles.push_back(clippedTriangle);
            }
        }
    }

    auto end1 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time for matrix operations: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count()
        << " ms" << std::endl;

    return processedTriangles;
}

void Engine3D::clipAndDraw(std::vector<Triangle>& triangles)
{
    auto start2 = std::chrono::steady_clock::now();

    // Clear the screen
    Clear(Pixel(0,0,0));

    // Sort the triangles based on the mean z to draw them in order.
    sort(triangles.begin(), triangles.end(), [](Triangle &t1, Triangle &t2)
    {
        float z1 = (t1.p[0].z() + t1.p[1].z() + t1.p[2].z()) / 3.0f;
        float z2 = (t2.p[0].z() + t2.p[1].z() + t2.p[2].z()) / 3.0f;
        return z1 > z2;
    });

    // Clip the triangles agains the four planes (X+, X-. Y+, Y-)
    for (auto& triangle : triangles)    // This should be parallelized
    {
        // List of triangles to clip
        std::deque<Triangle> trianglesToClip;

        trianglesToClip.push_back(triangle);
        int32_t nbNewTriangles = 1;

        // Clip all the triangles in the queue against a plane and add all the new triangles to the queue.
        // When all the triangles of the queue have been processed, go to the new plane.
        for (int32_t planeId = 0; planeId < 4; planeId++)
        {
            while (nbNewTriangles > 0)
            {
                Triangle test = trianglesToClip.front();
                trianglesToClip.pop_front();
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
                    trianglesToClip.push_back(clippedTriangle);
            }
            // All the triangles of the queue have been clipped agains the plane. Update the number of triangles to process for the new plane
            nbNewTriangles = trianglesToClip.size();
        }

        // Draw the triangles on the screen
        for (auto& triangleToRaster : trianglesToClip)
        {
            FillTriangle(triangleToRaster.p[0].x(), triangleToRaster.p[0].y(), 
                            triangleToRaster.p[1].x(), triangleToRaster.p[1].y(),  
                            triangleToRaster.p[2].x(), triangleToRaster.p[2].y(), Pixel(225, 25, 25));

            DrawTriangle(triangleToRaster.p[0].x(), triangleToRaster.p[0].y(), 
                            triangleToRaster.p[1].x(), triangleToRaster.p[1].y(),  
                            triangleToRaster.p[2].x(), triangleToRaster.p[2].y(), Pixel(255, 255, 255));

            std::stringstream textOverlay;
            textOverlay << "Rendering (ms) " << elapsedTime_ * 1000;
            DrawString(200, 0, textOverlay.str());
        }
    }

    auto end2 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time for geometric operations: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count()
                << " ms" << std::endl;
}


void Engine3D::updateCamera(void)
{
    if (GetKey(Key::UP).bHeld)
        vCameraPos_.y() += 100.0f * elapsedTime_;	    // Move the camera up

    if (GetKey(Key::DOWN).bHeld)
        vCameraPos_.y() -= 100.0f * elapsedTime_;	    // Move the camera down

    // Maybe add left right movement ???

    Vector3f vForwardMove = vLookDir_ * 100.0f * elapsedTime_;  // Forward movement if a key was pressed

    if (GetKey(Key::W).bHeld)
        vCameraPos_ += vForwardMove;    // Move the camera forward

    if (GetKey(Key::S).bHeld)
        vCameraPos_ -= vForwardMove;    // Move the camera backward

    if (GetKey(Key::A).bHeld)
        yaw_ += 1.0f * elapsedTime_;    // Rotate the camera left

    if (GetKey(Key::D).bHeld)
        yaw_ -= 1.0f * elapsedTime_;    // Rotate the camera right

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
