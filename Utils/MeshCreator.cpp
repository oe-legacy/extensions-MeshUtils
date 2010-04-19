// Mesh Creator.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Utils/MeshCreator.h>

#include <Geometry/Mesh.h>
#include <Geometry/Material.h>
#include <Geometry/GeometrySet.h>
#include <Resources/DataBlock.h>

#include <Logging/Logger.h>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {

        MeshPtr CreatePlane(float size, unsigned int detail, Vector<3, float> color){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;
            
            unsigned int points = d * d;
            
            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            Vector<3, float> normal = Vector<3, float>(0, 1, 0);
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, halfSize, j * unit - halfSize);
                    
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, normal);
                    colors->SetElement(index, color);
                }
            }

            unsigned int quads = detail * detail;
            unsigned int* i = new unsigned int[6 * quads];
            IndicesPtr indices = IndicesPtr(new Indices(6 * quads, i));

            unsigned int index = 0;
            for (unsigned int m = 0; m < detail; ++m){
                for (unsigned int n = 0; n < detail; ++n){
                    // Index the (i, j)'th quad of 2 triangles
                    i[index++] = m + n * d;
                    i[index++] = m + (n+1) * d;
                    i[index++] = m+1 + n * d;

                    i[index++] = m+1 + n * d;
                    i[index++] = m + (n+1) * d;
                    i[index++] = m+1 + (n+1) * d;
                }
            }
            
            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }
        
        MeshPtr CreateCube(float size, unsigned int detail, Vector<3, float> color){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;

            enum sides {TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5};
            
            unsigned int bottomOffset = d * d;
            unsigned int leftOffset = 2 * d * d;
            unsigned int rightOffset = 3 * d * d;
            unsigned int frontOffset = 4 * d * d;
            unsigned int backOffset = 5 * d * d;
            
            unsigned int points = 6 * d * d;

            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            // Top side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 1, 0));
                }
            }
            // Bottom geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = bottomOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, -halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, -1, 0));
                }
            }
            // Front geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = frontOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, j * unit - halfSize, halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 0, 1));
                }
            }
            // Back geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = backOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, j * unit - halfSize, -halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 0, -1));
                }
            }
            // Right side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = rightOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(halfSize, i * unit - halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(1, 0, 0));
                }
            }
            // Left side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = leftOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(-halfSize, i * unit - halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(-1, 0, 0));
                }
            }

            // colors
            for (unsigned int i = 0; i < colors->GetSize(); ++ i)
                colors->SetElement(i, color);

            IndicesPtr indices = IndicesPtr(new Indices(36 * detail * detail));
            unsigned int* i = indices->GetData();

            // Top side indices
            unsigned int index = 0;
            for (unsigned int k = 0; k < 6; ++k){
                unsigned int offset = k * d * d;
                for (unsigned int m = 0; m < detail; ++m){
                    for (unsigned int n = 0; n < detail; ++n){
                        // Index the (i, j)'th quad of 2 triangles
                        if (k == LEFT || k == TOP || k == BACK){
                            i[index++] = offset + m + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            i[index++] = offset + m+1 + n * d;
                            
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            i[index++] = offset + m+1 + (n+1) * d;
                        }else{
                            i[index++] = offset + m + n * d;
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m+1 + (n+1) * d;
                            i[index++] = offset + m + (n+1) * d;
                        }
                    }
                }
            }

            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }

        MeshPtr CreateSphere(float radius, unsigned int detail, Vector<3, float> color, bool inverted){
            MeshPtr mesh = CreateCube(radius, detail, color);
            DataBlock<3, float>* vertices = (DataBlock<3, float>*) mesh->GetGeometrySet()->GetVertices().get();
            DataBlock<3, float>* normals = (DataBlock<3, float>*) mesh->GetGeometrySet()->GetNormals().get();
            
            float inv = inverted ? -1 : 1;

            for (unsigned int i = 0; i < vertices->GetSize(); ++i){
                Vector<3, float> vert = vertices->GetElement(i);
                vert.Normalize();
                normals->SetElement(i, inv * vert);
                vertices->SetElement(i, inv * vert * radius);
            }
            
            return mesh;
        }        

    }
}
