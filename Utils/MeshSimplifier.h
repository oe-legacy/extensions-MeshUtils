// Mesh Simplifier.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#ifndef _OE_MESH_SIMPLFIER_H_
#define _OE_MESH_SIMPLFIER_H_

#include <boost/shared_ptr.hpp>

namespace OpenEngine {
    namespace Geometry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
    }
    namespace Utils {

        class MeshSimplifier {
        public:
            static Geometry::MeshPtr Simplify(Geometry::MeshPtr mesh);
        };

    }
}

#endif
