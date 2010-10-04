// Mesh Transformer.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Utils/MeshTransformer.h>
#include <Geometry/GeometrySet.h>
#include <Geometry/Mesh.h>
#include <Resources/IDataBlock.h>

#include <Logging/Logger.h>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {
        namespace MeshTransformer {

            MeshPtr Translate(MeshPtr mesh, Vector<3, float> move) {
                GeometrySetPtr newGeom = Translate(mesh->GetGeometrySet(), move);
                
                return MeshPtr(new Mesh(mesh->GetIndices(),
                                        mesh->GetType(),
                                        newGeom,
                                        mesh->GetMaterial(),
                                        mesh->GetIndexOffset(),
                                        mesh->GetDrawingRange()));
            }
            
            GeometrySetPtr Translate(GeometrySetPtr geom, Vector<3, float> move) {
                GeometrySetPtr clone = geom->Clone();

                IDataBlockPtr verts = clone->GetVertices();
                *verts += move;

                return clone;
            }

            MeshPtr Rotate(MeshPtr mesh, Quaternion<float> rotate) {
                GeometrySetPtr newGeom = Rotate(mesh->GetGeometrySet(), rotate);
                
                return MeshPtr(new Mesh(mesh->GetIndices(),
                                        mesh->GetType(),
                                        newGeom,
                                        mesh->GetMaterial(),
                                        mesh->GetIndexOffset(),
                                        mesh->GetDrawingRange()));
            }

            GeometrySetPtr Rotate(GeometrySetPtr geom, Quaternion<float> rotate) {
                GeometrySetPtr clone = geom->Clone();

                IDataBlockPtr verts = clone->GetVertices();
                if (verts != NULL)
                    for (unsigned int i = 0; i < verts->GetSize(); ++i){
                        Vector<3, float> elem;
                        verts->GetElement(i, elem);
                        elem = rotate.RotateVector(elem);
                        verts->SetElement(i, elem);
                    }

                IDataBlockPtr norms = clone->GetNormals();
                if (norms != NULL)
                    for (unsigned int i = 0; i < norms->GetSize(); ++i){
                        Vector<3, float> elem;
                        norms->GetElement(i, elem);
                        elem = rotate.RotateVector(elem);
                        norms->SetElement(i, elem);
                    }

                return clone;
            }

        }
    }
}
