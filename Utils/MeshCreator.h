// Mesh Creator.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#ifndef _OE_MESH_CREATOR_H_
#define _OE_MESH_CREATOR_H_

#include <Math/Vector.h>
#include <boost/shared_ptr.hpp>

using OpenEngine::Math::Vector;

namespace OpenEngine {
    namespace Geometry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
    }
    namespace Utils {

        namespace MeshCreator {

            Geometry::MeshPtr CreatePlane(float size, 
                                          unsigned int detail = 1, 
                                          Vector<3, float> color = Vector<3, float>());
            Geometry::MeshPtr CreateCube(float size, 
                                         unsigned int detail = 1, 
                                         Vector<3, float> color = Vector<3, float>());
            Geometry::MeshPtr CreateSphere(float radius, 
                                           unsigned int detail = 1, 
                                           Vector<3, float> color = Vector<3, float>(), 
                                           bool inverted = false);

        }

    }
}

#endif
