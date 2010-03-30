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
#include <Geometry/Material.h>
#include <Geometry/GeometrySet.h>
#include <Resources/DataBlock.h>
#include <Math/Matrix.h>
#include <Logging/Logger.h>

#include <algorithm>
#include <vector>
#include <limits.h>
#include <cstring>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {


        struct VertexPair {
            unsigned int v1;
            unsigned int v2;
            VertexAttr p;
            double error;
            
            bool operator<(const VertexPair& p){
                return this->error < p.error;
            }
            
            VertexPair(unsigned int vert1, unsigned int vert2, VertexAttr a, double err)
                : v1(vert1), v2(vert2), p(a), error(err) {}
        };

        MeshPtr CreatePlane(float size, Vector<3, float> color, unsigned int detail){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;
            
            unsigned int points = d * d;
            
            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            Vector<3, float> normal = Vector<3, float>(0, 1, 0);
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, halfSize, j * unit - halfSize);
                    
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, normal);
                    colors->SetElement(index, color);
                }
            }

            unsigned int quads = detail * detail;
            unsigned int* i = new unsigned int[6 * quads];
            IndicesPtr indices = IndicesPtr(new Indices(6 * quads, i));

            unsigned int index = 0;
            for (unsigned int m = 0; m < detail; ++m){
                for (unsigned int n = 0; n < detail; ++n){
                    // Index the (i, j)'th quad of 2 triangles
                    i[index++] = m + n * d;
                    i[index++] = m + (n+1) * d;
                    i[index++] = m+1 + n * d;

                    i[index++] = m+1 + n * d;
                    i[index++] = m + (n+1) * d;
                    i[index++] = m+1 + (n+1) * d;
                }
            }
            
            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }
        
        MeshPtr CreateCube(float size, Vector<3, float> color, unsigned int detail){
            unsigned int d = detail + 1;
            float halfSize = size / 2;
            float unit = size / detail;

            enum sides {TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5};
            
            unsigned int topOffset = 0;
            unsigned int bottomOffset = d * d;
            unsigned int leftOffset = 2 * d * d;
            unsigned int rightOffset = 3 * d * d;
            unsigned int frontOffset = 4 * d * d;
            unsigned int backOffset = 5 * d * d;
            
            unsigned int points = 6 * d * d;

            Float3DataBlockPtr vertices = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr normals = Float3DataBlockPtr(new DataBlock<3, float>(points));
            Float3DataBlockPtr colors = Float3DataBlockPtr(new DataBlock<3, float>(points));
            GeometrySetPtr geom = GeometrySetPtr(new GeometrySet(vertices, normals, Resources::IDataBlockList(), colors));

            // Top side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 1, 0));
                }
            }
            // Bottom geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = bottomOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, -halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, -1, 0));
                }
            }
            // Front geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = frontOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, j * unit - halfSize, halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 0, 1));
                }
            }
            // Back geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = backOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(i * unit - halfSize, j * unit - halfSize, -halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(0, 0, -1));
                }
            }
            // Right side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = rightOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(halfSize, i * unit - halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(1, 0, 0));
                }
            }
            // Left side geometry
            for (unsigned int i = 0; i < d; ++i){
                for (unsigned int j = 0; j < d; ++j){
                    unsigned int index = leftOffset + i + j * d;
                    Vector<3, float> vertex = Vector<3, float>(-halfSize, i * unit - halfSize, j * unit - halfSize);
                    vertices->SetElement(index, vertex);
                    normals->SetElement(index, Vector<3, float>(-1, 0, 0));
                }
            }

            // colors
            for (unsigned int i = 0; i < colors->GetSize(); ++ i)
                colors->SetElement(i, color);

            IndicesPtr indices = IndicesPtr(new Indices(36 * detail * detail));
            unsigned int* i = indices->GetData();

            // Top side indices
            unsigned int index = 0;
            for (unsigned int k = 0; k < 6; ++k){
                unsigned int offset = k * d * d;
                for (unsigned int m = 0; m < detail; ++m){
                    for (unsigned int n = 0; n < detail; ++n){
                        // Index the (i, j)'th quad of 2 triangles
                        if (k == LEFT || k == TOP || k == BACK){
                            i[index++] = offset + m + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            i[index++] = offset + m+1 + n * d;
                            
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            i[index++] = offset + m+1 + (n+1) * d;
                        }else{
                            i[index++] = offset + m + n * d;
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m + (n+1) * d;
                            
                            i[index++] = offset + m+1 + n * d;
                            i[index++] = offset + m+1 + (n+1) * d;
                            i[index++] = offset + m + (n+1) * d;
                        }
                    }
                }
            }

            return MeshPtr(new Mesh(indices, TRIANGLES, geom, MaterialPtr(new Material())));
        }
        
        MeshPtr CreateSphere(float radius, Vector<3, float> color, unsigned int detail){
            MeshPtr mesh = CreateCube(radius, color, detail);
            DataBlock<3, float>* vertices = (DataBlock<3, float>*) mesh->GetGeometrySet()->GetVertices().get();
            DataBlock<3, float>* normals = (DataBlock<3, float>*) mesh->GetGeometrySet()->GetNormals().get();
            
            for (unsigned int i = 0; i < vertices->GetSize(); ++i){
                Vector<3, float> vert = vertices->GetElement(i);
                vert.Normalize();
                normals->SetElement(i, vert);
                vertices->SetElement(i, vert * radius);
            }
            
            return mesh;
        }

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
            // Check to see that the properties of the mesh are
            // supported.
            if (mesh->GetType() != TRIANGLES)
                return mesh;
            if (mesh->GetGeometrySet()->GetVertices()->GetType() != FLOAT)
                return mesh;
            if (mesh->GetGeometrySet()->GetVertices()->GetDimension() != 3)
                return mesh;

            // There is no such thing as negative distances
            if (edgeMargin < 0)
                edgeMargin = 0;


            GeometrySetPtr geom = mesh->GetGeometrySet();
            DataBlock<3, float>* vertices = (DataBlock<3, float>*) geom->GetVertices().get();
            unsigned int points = vertices->GetSize();
            IndicesPtr indices = mesh->GetIndices();
            unsigned int triangles = indices->GetSize() / 3;
            unsigned int collapses = points - (points * reduction) / 100.0f;

            logger.info << "Points " << points << logger.end;
            logger.info << "Indices " << indices->GetSize() << logger.end;
            logger.info << "Triangles " << triangles << logger.end;
            logger.info << "Collapses to perform " << collapses << logger.end;

            // Create and calculate quadrics
            Matrix<4, 4, double> quadrics[points];
            
            // Initialize all quadrics to the 0 matrix.
            for (unsigned int i = 0; i < points; ++i)
                quadrics[i] = Matrix<4, 4, double>(0.0f);

            // Calculate quadric errors for each vertex
            for (unsigned int i = 0; i < indices->GetSize(); i += 3){
                unsigned int i1 = (*indices)[i][0];
                Vector<3, double> v1 = (*vertices)[i1].ToDouble();
                unsigned int i2 = (*indices)[i+1][0];
                Vector<3, double> v2 = (*vertices)[i2].ToDouble();
                unsigned int i3 = (*indices)[i+2][0];
                Vector<3, double> v3 = (*vertices)[i3].ToDouble();
                
                Vector<3, double> span1 = v2 - v1;
                Vector<3, double> span2 = v3 - v1;
                
                Vector<3, double> normal = (span1 % span2);
                normal.Normalize();
                //logger.info << "Normal " << normal << logger.end;

                float a = normal[0];
                float b = normal[1];
                float c = normal[2];
                float d = 0 - normal*v1;

                Matrix<4, 4, double> quad = 
                    Matrix<4, 4, double>(a*a, a*b, a*c, a*d,
                                        b*a, b*b, b*c, b*d,
                                        c*a, c*b, c*c, c*d,
                                        d*a, d*b, d*c, d*d);

                quadrics[i1] = quadrics[i1] + quad;
                quadrics[i2] = quadrics[i2] + quad;
                quadrics[i3] = quadrics[i3] + quad;
            }

            /**
             * For each point create a list of it's neighbours ...
             */
            vector<list<unsigned int> > neighbours = vector<list<unsigned int> >(points);
            for (unsigned int i = 0; i < indices->GetSize(); i += 3){
                unsigned int i1 = (*indices)[i][0];
                unsigned int i2 = (*indices)[i+1][0];
                unsigned int i3 = (*indices)[i+2][0];
                
                neighbours[i1].push_back(i2);
                neighbours[i1].push_back(i3);
                neighbours[i2].push_back(i1);
                neighbours[i2].push_back(i3);
                neighbours[i3].push_back(i1);
                neighbours[i3].push_back(i2);
            }
            // ... and then sort it.
            for (unsigned int i = 0; i < points; ++i){
                neighbours[i].sort();
            }

            /**
             * Create the collapsable vector of vertice attributes.
             */
            vector<VertexAttr> verticeAttr;
            for (unsigned int i = 0; i < points; ++i){
                verticeAttr.push_back(CreateVertexAttr(geom, i));
            }

            vector<list<VertexPair*> > vertexPairs = vector<list<VertexPair*> >(points);

            /**
             * Create a list of pairs that can be collapsed.
             */
            list<VertexPair*> collapsables;
            for (unsigned int i = 0; i < points; ++i){
                Vector<3, float> v = (*vertices)[i];
                Vector<4, float> v1 = Vector<4, float>(v[0], v[1], v[2], 1);

                list<unsigned int>::iterator nStart = neighbours[i].begin();
                list<unsigned int>::iterator nEnd = neighbours[i].end();

                for (unsigned int j = i; j < points; ++j){
                    Vector<3, float> v = (*vertices)[j];
                    Vector<4, float> v2 = Vector<4, float>(v[0], v[1], v[2], 1);
                
                    // If either the vertices are neighbours or the
                    // distance between them is low enough the pair
                    // can be considered for contraction.
                    if (find(nStart, nEnd, j) != nEnd || 
                        (v1 - v2).GetLength() < edgeMargin) {
                        
                        Matrix<4, 4, double> quad = quadrics[i] + quadrics[j];
                        float bias1 = 0.5f;
                        float bias2 = 0.5f;
                        VertexAttr attr = verticeAttr[i] * bias1 + verticeAttr[j] * bias2;
                        
                        double error = (quad * attr.vec.ToDouble()) * attr.vec.ToDouble();
                        
                        VertexPair* pair = new VertexPair(i, j, attr, error);

                        collapsables.push_back(pair);

                        vertexPairs[i].push_back(pair);
                        vertexPairs[j].push_back(pair);
                    }
                }
            }
            neighbours.clear();

            /**
             * Sort the list in increasing order, so we can start
             * collapsing the least important pairs first.
             */
            collapsables.sort();

            /**
             * Collapse the pairs. After each collapse the new vertex
             * is inserted into the pair list and the indices is
             * updated. Also a flag is set indicating the vertex that
             * was removed.
             */
            unsigned int collapsed[points];
            for (unsigned int i = 0; i < points; ++i)
                collapsed[i] = UINT_MAX;

            IndicesPtr newIs = IndicesPtr(new Indices(indices->GetSize()));
            memcpy(newIs->GetData(), indices->GetData(), indices->GetSize());

            while (collapses > 0 && collapsables.size() > 0){
                VertexPair* pair = collapsables.front();
                unsigned int updatedIndex = pair->v1;
                unsigned int deletedIndex = pair->v2;

                // Update geometry
                verticeAttr[updatedIndex] = pair->p;

                // Update neighbours neighbours

                // Update errors

                // Update indices in the new indices, the pairs and the collapsables.

                collapsed[deletedIndex] = updatedIndex;
                collapsables.pop_front();
                --collapses;
            }

            /**
             * Create a new mesh from the collapsed vertices.
             */
            

            /**
             * Remove degenerate triangles from the indices.
             */
            

            return mesh;
        }

        /**
         * Creates a VertexAttr from the i'th attributes in the GeometrySet.
         */
        VertexAttr CreateVertexAttr(GeometrySetPtr geom, unsigned int i){
            Vector<4, float> v = Vector<4, float>((float*)geom->GetVertices()->GetVoidElement(i));
            Vector<3, float> n = Vector<3, float>((float*)geom->GetNormals()->GetVoidElement(i));
            Vector<4, float> c = Vector<4, float>((float*)geom->GetColors()->GetVoidElement(i));
            list<Vector<3, float> > tc;
            list<IDataBlockPtr>::iterator itr = geom->GetTexCoords().begin();
            while(!geom->GetTexCoords().empty() && itr != geom->GetTexCoords().end()){
                tc.push_back(Vector<3, float>((float*)(*itr)->GetVoidElement(i)));
                ++itr;
            }
            return VertexAttr(v, n, c, tc);
        }

        /**
         * Eliminates degenerate triangles from a list of indices.
         *
         * @param i The list of indices to be reduced.
         *
         * @return the new list of indices.
         */
        IndicesPtr RemoveDegenerates(IndicesPtr i){
            unsigned int* tmp = new unsigned int[i->GetSize()];
            unsigned int offset = 0;
            unsigned int k = 0;            
            for (k = 0; k < i->GetSize(); k += 3) {
                unsigned int i1 = i->GetData()[k];
                unsigned int i2 = i->GetData()[k+1];
                unsigned int i3 = i->GetData()[k+2];
                if (i1 != i2 && i1 != i3 && i2 != i3){
                    tmp[offset++] = i1;
                    tmp[offset++] = i2;
                    tmp[offset++] = i3;
                }
            }
            
            unsigned int* indices = new unsigned int[k];
            memcpy(indices, tmp, offset * sizeof(unsigned int));
            return IndicesPtr(new Indices(k, indices));
        }
        
    }
}
