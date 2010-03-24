// Mesh Utils.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Utils/MeshUtils.h>
#include <Geometry/Mesh.h>
#include <Geometry/Material.h>
#include <Geometry/GeometrySet.h>
#include <Resources/DataBlock.h>
#include <Logging/Logger.h>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {

        MeshPtr CreatePlane(float size, Vector<3, float> color, unsigned int detail){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;
            
            unsigned int points = d * d;
            
            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            unsigned int quads = detail * detail;
            unsigned int* i = new unsigned int[6 * quads];
            DataIndicesPtr indices = DataIndicesPtr(new DataIndices(6 * quads, i));

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
            unsigned int index = 0;
            for (unsigned int m = 0; m < detail; ++m){
                for (unsigned int n = 0; n < detail; ++n){
                    // Index the (i, j)'th quad of 2 triangles
                    i[index] = m + n * d;
                    ++index;
                    i[index] = m+1 + n * d;
                    ++index;
                    i[index] = m + (n+1) * d;
                    ++index;

                    i[index] = m+1 + n * d;
                    ++index;
                    i[index] = m + (n+1) * d;
                    ++index;
                    i[index] = m+1 + (n+1) * d;
                    ++index;
                }
            }
            
            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }
        
        MeshPtr CreateCube(float size, unsigned int detail, Vector<3, float> color){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;
            
            unsigned int topOffset = 0;
            unsigned int bottomOffset = d * d;
            unsigned int leftOffset = 2 * d * d;
            unsigned int rightOffset = 3 * d * d;
            unsigned int frontOffset = 4 * d * d;
            unsigned int backOffset = 5 * d * d;
            
            unsigned int points = 6 * d * d;
            float* v = new float[3 * points];
            float* n = new float[3 * points];
            float* c = new float[3 * points];
            
            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points, v));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points, n));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points, c));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            unsigned int* i = new unsigned int[6 * points];
            DataIndicesPtr indices = DataIndicesPtr(new DataIndices(6 * points, i));

            // Top side geometry
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
            // Top side indices
            unsigned int index = 0;
            for (unsigned int m = 0; m < detail; ++m){
                for (unsigned int n = 0; n < detail; ++n){
                    // Index the (i, j)'th quad of 2 triangles
                    i[index] = topOffset + m + n * d;
                    ++index;
                    i[index] = topOffset + m+1 + n * d;
                    ++index;
                    i[index] = topOffset + m + (n+1) * d;
                    ++index;

                    i[index] = topOffset + m+1 + n * d;
                    ++index;
                    i[index] = topOffset + m + (n+1) * d;
                    ++index;
                    i[index] = topOffset + m+1 + (n+1) * d;
                    ++index;
                }
            }
            
            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }
        
        MeshPtr CreateSphere(float radius, unsigned int detail, Vector<3, float> color){
            return MeshPtr();
        }

        /**
         * This method takes a mesh and simplifies it, meaning making
         * it less detailed. Only works if all datablocks are of the
         * type float.
         *
         * @param mesh The mesh to simplify.
         * @param onlyEdges Decides wether to only contract linked
         * edges or also contract non linked edges.
         * @param reduction How many percent should the number of
         * vertices in the model be be reduced by.
         *
         * @return The simplifies mesh.
         */
        MeshPtr Simplify(MeshPtr mesh, float edgeMargin, char reduction) {
#ifdef OE_SAFE
            if (mesh->GetPrimitive() != TRIANGLES)
                throw Exception("Unsupported geometry primitive.");
#endif
            if (edgeMargin < 0)
                edgeMargin = 0;
            

            return mesh;
        }

    }
}
