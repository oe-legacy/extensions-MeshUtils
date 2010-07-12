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
            double error;
            
            bool operator<(const VertexPair& p){
                return this->error < p.error;
            }
            
            VertexPair(unsigned int vert1, unsigned int vert2, double err)
                : v1(vert1), v2(vert2), error(err) {}
        };

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
            if (mesh->GetGeometrySet()->GetVertices()->GetType() != Types::FLOAT)
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
            unsigned int maxContractions = points - (points * reduction) / 100.0f;

            logger.info << "Points " << points << logger.end;
            logger.info << "Indices " << indices->GetSize() << logger.end;
            logger.info << "Contractions to perform " << maxContractions << logger.end;

            // Create and calculate quadrics
            vector<Matrix<4, 4, double> > quadrics(points);
            
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
            vector<VertexAttr> vertexAttrs;
            for (unsigned int i = 0; i < points; ++i){
                vertexAttrs.push_back(CreateVertexAttr(geom, i));
            }

            /**
             * A reference to the pairs that the i'th vertex is 
             */
            vector<list<VertexPair*> > vertexPairs = vector<list<VertexPair*> >(points);
            
            /**
             * Create a list of pairs that can be contracted.
             */
            list<VertexPair*> contractables;
            for (unsigned int i = 0; i < points; ++i){
                Vector<3, float> v = (*vertices)[i];
                Vector<4, double> v1 = Vector<4, double>(v[0], v[1], v[2], 1);
                
                list<unsigned int>::iterator nStart = neighbours[i].begin();
                list<unsigned int>::iterator nEnd = neighbours[i].end();
                
                for (unsigned int j = i; j < points; ++j){
                    Vector<3, float> v = (*vertices)[j];
                    Vector<4, double> v2 = Vector<4, double>(v[0], v[1], v[2], 1);
                    
                    // If either the vertices are neighbours or the
                    // distance between them is low enough the pair
                    // can be considered for contraction.
                    if (find(nStart, nEnd, j) != nEnd || 
                        (v1 - v2).GetLength() < edgeMargin) {

                        // The error term can be mathematically
                        // reduced to a lookup into the inverse error
                        // matrix.
                        
                        // double error = quad.GetInverse()(3, 3);
                        
                        // This of course has to handle the special
                        // case where the matrix isn't invertible. The
                        // error then is probably 0 and a good new
                        // vector is (v1 + v2 / 2).
                        
                        Matrix<4, 4, double> quad = quadrics[i] + quadrics[j];
                        Vector<4, double> vec = ((v1 + v2) / 2.0).ToDouble();
                        
                        double error = (quad * vec) * vec;
                        
                        VertexPair* pair = new VertexPair(i, j, error);
                        
                        contractables.push_back(pair);
                        
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
            contractables.sort();
            
            /**
             * Collapse the pairs. After each collapse the new vertex
             * is inserted into the pair list and the indices is
             * updated. Also the vertex removed must point to it's new
             * index.
             */
            unsigned int contractedTo[points];
            for (unsigned int i = 0; i < points; ++i)
                contractedTo[i] = UINT_MAX;
            
            IndicesPtr newIs = IndicesPtr(new Indices(indices->GetSize()));
            memcpy(newIs->GetData(), indices->GetData(), indices->GetSize());
            
            unsigned int contractions = 0;
            while (contractions < maxContractions && 0 < contractables.size()){
                //logger.info << "Contractions performed " << contractions << logger.end;
                VertexPair* pair = contractables.front();
                //logger.info << pair->v1 << " and " << pair->v2 << logger.end;
                unsigned int i = GetContractedIndex(pair->v1, contractedTo);
                unsigned int j = GetContractedIndex(pair->v2, contractedTo);
                //logger.info << pair->v1 << " has been moved to " << i << logger.end;
                //logger.info << pair->v2 << " has been moved to " << j << logger.end;
                
                if (i != j){

                    // Sanity check
                    if (vertexAttrs[i].live == false ||
                        vertexAttrs[j].live == false)
                        logger.info << "WTF! A dead vertex made it through." << logger.end;
                
                    Matrix<4, 4, double> quadric = quadrics[i] + quadrics[j];
                    
                    // New vertex attr
                    Vector<4, double> vec;
                    if (quadric.GetDeterminant() != 0.0)
                        vec = quadric.GetInverse() * Vector<4, double>(0, 0, 0, 1);
                    else{
                        Vector<3, float> v = (*vertices)[i] + (*vertices)[j];
                        vec = Vector<4, double>(v[0], v[1], v[2], 0);
                    }
                    
                    // Contract the vertex attributes based on a bias
                    // calculated from vec.
                    float bias1 = (Vector<3, float>(vec[0], vec[1], vec[2]) - (*vertices)[i]).GetLengthSquared();
                    float bias2 = (Vector<3, float>(vec[0], vec[1], vec[2]) - (*vertices)[j]).GetLengthSquared();
                    vertexAttrs[i] *= bias1;
                    vertexAttrs[j] *= bias2;
                    vertexAttrs[i] += vertexAttrs[j];
                    vertexAttrs[i].vec = Vector<4, float>(vec[0], vec[1], vec[2], vec[3]);
                    
                    // Update the quadric errors for the new Vertex
                    quadrics[i] = quadric;

                    // Update the contracted to reference and number of contractions.
                    contractedTo[j] = i;
                    vertexAttrs[j].live = false;
                    contractions++;
                }

                contractables.pop_front();
            }
                
            /**
             * Create a new mesh from the collapsed vertices.
             */
            GeometrySetPtr newGeom = CreateGeometry(geom, vertexAttrs, contractions);
            
            /**
             * Update the indices based on information in
             * contractedTo and then remove degenerate triangles.
             */
            logger.info << "Original indices " << indices->ToString() << logger.end;
            IndicesPtr newIndices = UpdateIndices(indices, contractedTo);
            logger.info << "New indices " << newIndices->ToString() << logger.end;
            newIndices = RemoveDegenerates(newIndices);
            logger.info << "New indices " << newIndices << logger.end;

            return mesh;
            /*
            return MeshPtr(new Mesh(newIndices, 
                                    TRIANGLES, 
                                    geom, 
                                    mesh->GetMaterial()));
            */
        }

        /**
         * Creates a VertexAttr from the i'th attributes in the GeometrySet.
         */
        VertexAttr CreateVertexAttr(GeometrySetPtr geom, unsigned int i){
            Vector<4, float> v;// = Vector<4, float>((float*)geom->GetVertices()->GetVoidElement(i));
            geom->GetVertices()->GetElement(i, v);
            Vector<3, float> n;// = Vector<3, float>((float*)geom->GetNormals()->GetVoidElement(i));
            geom->GetNormals()->GetElement(i, n);
            Vector<4, float> c;// = Vector<4, float>((float*)geom->GetColors()->GetVoidElement(i));
            geom->GetColors()->GetElement(i, c);
            list<Vector<3, float> > tc;
            list<IDataBlockPtr>::iterator itr = geom->GetTexCoords().begin();
            while(!geom->GetTexCoords().empty() && itr != geom->GetTexCoords().end()){
                Vector<3, float> t;
                (*itr)->GetElement(i, t);
                //tc.push_back(Vector<3, float>((float*)(*itr)->GetVoidElement(i)));
                tc.push_back(t);
                ++itr;
            }
            return VertexAttr(v, n, c, tc);
        }

        /**
         * If the vertex at index i, has been contracted into another
         * index, then contractedTo is recursively searched until the
         * final index is found.
         *
         * @param i The original index.
         * @param contractedTo An array of indices that an index has
         * been contracted into.
         *
         * @return The index of the vertex after the previous contractions.
         */
        unsigned int GetContractedIndex(unsigned int i, unsigned int* contractedTo){
            unsigned int j = i;
            while (contractedTo[j] != UINT_MAX)
                j = contractedTo[j];
            
            return j;
        }
            
        /**
         * Updates the indices based on which new index the original
         * index was contracted to.
         */
        IndicesPtr UpdateIndices(IndicesPtr i, unsigned int* contractedTo){
            IndicesPtr ret = IndicesPtr(new Indices(i->GetSize()));
            for (unsigned int j = 0; j < i->GetSize(); ++j)
                ret->GetData()[j] = GetContractedIndex(i->GetData()[j], contractedTo);
            
            return ret;
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
            return IndicesPtr(new Indices(offset, indices));
        }
        
        GeometrySetPtr CreateGeometry(GeometrySetPtr original, vector<VertexAttr> attrs, unsigned int size){
            IDataBlockPtr vertices = IDataBlockPtr();
            if (original->GetVertices() != NULL){
                unsigned int size = original->GetVertices()->GetSize();
                if (original->GetVertices()->GetDimension() == 3){
                    Float3DataBlockPtr vs = Float3DataBlockPtr(new DataBlock<3, float>(size));
                    unsigned int o = 0;
                    for (unsigned int i = 0; i < size; ++i){
                        if (attrs[i+o].live){
                            Vector<4, float> v = attrs[i].vec;
                            Vector<3, float> vec = Vector<3, float>(v[0],v[1],v[2]);
                            vs->SetElement(i, vec);
                        }else
                            o += 1;
                    }
                    vertices = vs;
                    
                }else if (original->GetVertices()->GetDimension() == 4){
                    Float4DataBlockPtr vs = Float4DataBlockPtr(new DataBlock<4, float>(size));
                    unsigned int o = 0;
                    for (unsigned int i = 0; i < size; ++i){
                        if (attrs[i+o].live){
                            vs->SetElement(i, attrs[i].vec);
                        }else
                            o += 1;
                    }
                    vertices = vs;
                }
            }
            
            return GeometrySetPtr();
        }
        
    }
}
