// Mesh Simplifier.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Utils/MeshSimplifier.h>

#include <Geometry/Mesh.h>
#include <Geometry/GeometrySet.h>
#include <Geometry/Vertex.h>

#include <vector>

using namespace OpenEngine::Geometry;
using namespace std;

namespace OpenEngine {
    namespace Utils {
        
        MeshPtr MeshSimplifier::Simplify(MeshPtr mesh){
            unsigned int oldSize = mesh->GetGeometrySet()->GetSize();
            vector<Vertex<double> > vertices(oldSize);
            for (unsigned int i = 0; i < oldSize; ++i){
                vertices[i] = mesh->GetGeometrySet()->GetVertex<double>(i);
            }

            return mesh;
        }

    }
}
