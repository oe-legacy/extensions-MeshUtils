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

using namespace OpenEngine::Geometry;

namespace OpenEngine {
    namespace Utils {

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
