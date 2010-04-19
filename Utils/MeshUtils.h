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

#include <list>
#include <vector>

using namespace std;

using namespace OpenEngine::Math;

namespace OpenEngine {
    namespace Geometry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
        class GeometrySet;
        typedef boost::shared_ptr<GeometrySet> GeometrySetPtr;
    }
    namespace Resources {
        class Indices;
        typedef boost::shared_ptr<Indices> IndicesPtr;
    }
    namespace Utils {

        /**
         * Is this the basis for a point class?
         * Could be useful both here and in a triangle/face class.
         */
        struct VertexAttr {
            Vector<4, float> vec;
            Vector<3, float> norm;
            Vector<4, float> color;
            list<Vector<3, float> > texCoords;
            bool live;

            VertexAttr operator*=(const float s){
                vec *= s;
                norm *= s;
                color *= s;
                list<Vector<3, float> >::iterator itr = texCoords.begin();
                while(itr != texCoords.end()){
                    (*itr) *= s;
                    ++itr;
                }
                return *this;
            }

            VertexAttr operator+=(const VertexAttr a){
                vec += a.vec;
                norm += a.norm;
                color += a.color;
                list<Vector<3, float> >::iterator itr1 = texCoords.begin();
                list<Vector<3, float> >::const_iterator itr2 = a.texCoords.begin();
                while(itr1 != texCoords.end()){
                    (*itr1) += (*itr1);
                    ++itr1; ++itr2;
                }
                return *this;
            }
            
            VertexAttr(Vector<4, float> v, Vector<3, float> n, Vector<4, float> c, list<Vector<3, float> > t)
                : vec(v), norm(n), color(c), texCoords(t) {
                live = true;
            }
        };

        Geometry::MeshPtr Simplify(Geometry::MeshPtr mesh, float edgeMargin = 0, char reduction = 75);
        VertexAttr CreateVertexAttr(Geometry::GeometrySetPtr geom, unsigned int i);
        unsigned int GetContractedIndex(unsigned int i, unsigned int* contractedTo);
        Resources::IndicesPtr UpdateIndices(Resources::IndicesPtr i, unsigned int* contractedTo);
        Resources::IndicesPtr RemoveDegenerates(Resources::IndicesPtr i);
        Geometry::GeometrySetPtr CreateGeometry(Geometry::GeometrySetPtr original, vector<VertexAttr> attrs, unsigned int contractions);
    }
}

#endif
