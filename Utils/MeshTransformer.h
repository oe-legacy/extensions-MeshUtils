// Mesh Transformer.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#ifndef _OE_MESH_TRANSFORMER_H_
#define _OE_MESH_TRANSFORMER_H_

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

            Geometry::MeshPtr Translate(Geometry::MeshPtr mesh, Vector<3, float> move);

            /**
             * Remember to rotate the normals aswell when rotating.
             */

            /**
             * Remember to scale the normals aswell when scaling and
             * then normalizing them again.
             */
        }
    }
}

#endif
