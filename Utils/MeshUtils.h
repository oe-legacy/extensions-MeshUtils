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

namespace OpenEngine {
    namespace Geoemtry {
        class Mesh;
        typedef boost::shared_ptr<Mesh> MeshPtr;
    }
    namespace Utils {
        
        MeshPtr Simplify(MeshPtr mesh, float edgeMargin = 0, char reduction = 75){
            
        }

    }
}

#endif
