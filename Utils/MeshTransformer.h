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
#include <Math/Quaternion.h>
#include <boost/shared_ptr.hpp>

using namespace OpenEngine::Math;

namespace OpenEngine {
    namespace Geometry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
        class GeometrySet;
        typedef boost::shared_ptr<GeometrySet> GeometrySetPtr;
    }
    namespace Utils {

        namespace MeshTransformer {

            /**
             * Build transformers that can modify an entire subtree.
             */

            Geometry::MeshPtr Translate(Geometry::MeshPtr geom, Vector<3, float> move);
            Geometry::GeometrySetPtr Translate(Geometry::GeometrySetPtr geom, Vector<3, float> move);

            /**
             * Remember to rotate the normals aswell when rotating.
             */
            Geometry::MeshPtr Rotate(Geometry::MeshPtr geom, Quaternion<float> rotate);
            Geometry::GeometrySetPtr Rotate(Geometry::GeometrySetPtr geom, Quaternion<float> rotate);

            /**
             * Remember to scale the normals aswell when scaling and
             * then normalizing them again. Needs to be done for non
             * uniform scaling.
             */
        }
    }
}

#endif
