// Mesh Utils.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#ifndef _OE_MESH_UTILS_H_
#define _OE_MESH_UTILS_H_

#include <boost/shared_ptr.hpp>
#include <Math/Vector.h>

using namespace OpenEngine::Math;

namespace OpenEngine {
    namespace Geometry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
    }
    namespace Utils {

        Geometry::MeshPtr CreatePlane(float size, Vector<3, float> color, unsigned int detail = 1);
        Geometry::MeshPtr CreateCube(float size, unsigned int detail, Vector<3, float> color);
        Geometry::MeshPtr CreateSphere(float radius, unsigned int detail, Vector<3, float> color);
        
        Geometry::MeshPtr Simplify(Geometry::MeshPtr mesh, float edgeMargin = 0, char reduction = 75);

    }
}

#endif
