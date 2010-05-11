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
#include <Resources/IDataBlock.h>

#include <Logging/Logger.h>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {
        namespace MeshTransformer {
            
            GeometrySetPtr Translate(GeometrySetPtr geom, Vector<3, float> move) {
                GeometrySetPtr clone = geom->Clone();

                IDataBlockPtr verts = clone->GetVertices();
                *verts += move;

                return clone;
            }

        }
    }
}
