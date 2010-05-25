// Mesh Creator.
// -------------------------------------------------------------------
// Copyright (C) 2010 OpenEngine.dk (See AUTHORS) 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Utils/MeshCreator.h>

#include <Geometry/Mesh.h>
#include <Geometry/Material.h>
#include <Geometry/GeometrySet.h>
#include <Resources/DataBlock.h>

#include <Logging/Logger.h>

#include<vector>

using namespace OpenEngine::Geometry;
using namespace OpenEngine::Resources;

namespace OpenEngine {
    namespace Utils {
        namespace MeshCreator {

            MeshPtr CreatePlane(float size, unsigned int detail,
                                Vector<3, float> color){
                unsigned int d = detail + 1;
                float halfSize = size / 2;
                float unit = size / detail;
            
                unsigned int points = d * d;
            
                Float3DataBlockPtr vertices =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                Float3DataBlockPtr normals =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                IDataBlockList texCoordList;            
                Float2DataBlockPtr texCoords = 
                    Float2DataBlockPtr(new DataBlock<2, float>(points));
                texCoordList.push_back(texCoords);

                Float3DataBlockPtr colors = 
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                GeometrySetPtr geom = 
                    GeometrySetPtr(new GeometrySet(vertices, 
                                                   normals, 
                                                   texCoordList,
                                                   colors));

                Vector<3, float> normal = Vector<3, float>(0, 1, 0);
                for (unsigned int i = 0; i < d; ++i){
                    for (unsigned int j = 0; j < d; ++j){
                        unsigned int index = i + j * d;
                        Vector<3, float> vertex(i * unit - halfSize,
                                                0, 
                                                j * unit - halfSize);
                    
                        vertices->SetElement(index, vertex);
                        normals->SetElement(index, normal);
                        Vector<2,float> texCoord( ((float)i)/detail,
                                                  ((float)j)/detail);
                        texCoords->SetElement(index, texCoord);
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
            
            return MeshPtr(new Mesh(indices, TRIANGLES, 
                                    geom, MaterialPtr(new Material())));
            }
        
            MeshPtr CreateCube(float size, unsigned int detail, Vector<3, float> color){
                unsigned int d = detail + 1;
                float halfSize = size / 2;
                float unit = size / detail;

                enum sides {TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5};
            
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

                return MeshPtr(new Mesh(indices, TRIANGLES, geom, 
                                        MaterialPtr(new Material())));
            }

            MeshPtr CreateSphere(float radius, unsigned int detail,
                                 Vector<3, float> color, bool inverted){
                MeshPtr mesh = CreateCube(radius, detail, color);
                DataBlock<3, float>* vertices = (DataBlock<3, float>*)
                    mesh->GetGeometrySet()->GetVertices().get();
                DataBlock<3, float>* normals = (DataBlock<3, float>*) 
                    mesh->GetGeometrySet()->GetNormals().get();
            
                float inv = inverted ? -1 : 1;

                for (unsigned int i = 0; i < vertices->GetSize(); ++i){
                    Vector<3, float> vert = vertices->GetElement(i);
                    vert.Normalize();
                    normals->SetElement(i, inv * vert);
                    vertices->SetElement(i, inv * vert * radius);
                }
            
                return mesh;
            }        

            class Face3 {
            public:
                int ids[3];
                Face3() {}
            };
            
            typedef Vector<3,float> Vector3;
            // from: http://www.google.com/codesearch/p?hl=en#IEMqheNvKao/trunk/guitest/b3d_lib/src/geodesic.cpp&q=Geodesic%20sphere%20lang:c%2B%2B&sa=N&cd=8&ct=rc&l=3
            MeshPtr create(int subdivisions, double radius,
                           int inverted, Vector<3,float> color) {
                // create a icosahedron first
                //double t1 = 2*PI/5;
                double t2 = PI/10;
                //double t3 = -3*PI/10;
                double t4 = PI/5;

                //double R = (side/2)/sin(t4);
                //double side = 1.0;
                double R = radius;
                double side = 2*R*sin(t4);
                double S = side;

                double H = cos(t4)*R;
                double Cx = R * cos(t2);
                double Cy = R * sin(t2);

                double H1 = sqrt(side*side - R*R);
                double H2 = sqrt((H+R)*(H+R) - H*H);
                double Z2 = (H2 - H1)/2;
                double Z1 = Z2 + H1;

                std::vector<Vector3> pts(12);
                pts[0] = Vector3(   0,  0, Z1);
                pts[1] = Vector3(   0,  R, Z2);
                pts[2] = Vector3(  Cx, Cy, Z2);
                pts[3] = Vector3( S/2, -H, Z2);
                pts[4] = Vector3(-S/2, -H, Z2);
                pts[5] = Vector3( -Cx, Cy, Z2);
                pts[6] = Vector3(   0, -R,-Z2);
                pts[7] = Vector3( -Cx,-Cy,-Z2);
                pts[8] = Vector3(-S/2,  H,-Z2);
                pts[9] = Vector3( S/2,  H,-Z2);
                pts[10] = Vector3( Cx,-Cy,-Z2);
                pts[11] = Vector3(  0,  0,-Z1);

                // now create the faces
                std::vector<Face3> face(20);
                face[0].ids[0] = 0;
                face[0].ids[1] = 1;
                face[0].ids[2] = 2;
                face[1].ids[0] = 0;
                face[1].ids[1] = 2;
                face[1].ids[2] = 3;
                face[2].ids[0] = 0;
                face[2].ids[1] = 3;
                face[2].ids[2] = 4;
                face[3].ids[0] = 0;
                face[3].ids[1] = 4;
                face[3].ids[2] = 5;
                face[4].ids[0] = 0;
                face[4].ids[1] = 5;
                face[4].ids[2] = 1;
                face[5].ids[0] = 1;
                face[5].ids[1] = 9;
                face[5].ids[2] = 2;
                face[6].ids[0] = 2;
                face[6].ids[1] = 9;
                face[6].ids[2] = 10;
                face[7].ids[0] = 2;
                face[7].ids[1] = 10;
                face[7].ids[2] = 3;
                face[8].ids[0] = 3;
                face[8].ids[1] = 10;
                face[8].ids[2] = 6;
                face[9].ids[0] = 3;
                face[9].ids[1] = 6;
                face[9].ids[2] = 4;
                face[10].ids[0] = 4;
                face[10].ids[1] = 6;
                face[10].ids[2] = 7;
                face[11].ids[0] = 4;
                face[11].ids[1] = 7;
                face[11].ids[2] = 5;
                face[12].ids[0] = 5;
                face[12].ids[1] = 7;
                face[12].ids[2] = 8;
                face[13].ids[0] = 5;
                face[13].ids[1] = 8;
                face[13].ids[2] = 1;
                face[14].ids[0] = 1;
                face[14].ids[1] = 8;
                face[14].ids[2] = 9;

                face[15].ids[0] = 9;
                face[15].ids[1] = 11;
                face[15].ids[2] = 10;
                face[16].ids[0] = 10;
                face[16].ids[1] = 11;
                face[16].ids[2] = 6;
                face[17].ids[0] = 6;
                face[17].ids[1] = 11;
                face[17].ids[2] = 7;
                face[18].ids[0] = 7;
                face[18].ids[1] = 11;
                face[18].ids[2] = 8;
                face[19].ids[0] = 8;
                face[19].ids[1] = 11;
                face[19].ids[2] = 9;

                // now we have our icosahedron

                // let's extend it to a geodesic sphere
                // subdivide each face
                std::vector<Face3> facesGS = face;
                std::vector<Vector3> nodesGS = pts;
                for (int subdiv=1; subdiv<=subdivisions; ++subdiv) {
                    // temp face list (previous subdivision)
                    std::vector<Face3> facesT = facesGS;    
                    
                    // temp node list
                    std::vector<Vector3> nodesT = nodesGS;
                    
                    facesGS.clear();
                    //nodesGS.clear();
                    for (size_t i=0; i<facesT.size(); ++i) {
                        Face3 &f = facesT[i];
                        
                        // create the new points
                        Vector3 newV1 = 
                            (nodesT[f.ids[0]] + nodesT[f.ids[1]])/2.0;
                        Vector3 newV2 =
                            (nodesT[f.ids[1]] + nodesT[f.ids[2]])/2.0;
                        Vector3 newV3 =
                            (nodesT[f.ids[2]] + nodesT[f.ids[0]])/2.0;
                        
                        // push the new points out to the radius of the sphere
                        //double len = newV1.length();
                          //newV1.x = newV1.x * R / len;
                          // newV1.y = newV1.y * R / len;
//                           newV1.z = newV1.z * R / len;
//                         len = newV2.length();
//                         newV2.x = newV2.x * R / len;
//                         newV2.y = newV2.y * R / len;
//                         newV2.z = newV2.z * R / len;
//                         len = newV3.length();
//                         newV3.x = newV3.x * R / len;
//                         newV3.y = newV3.y * R / len;
//                         newV3.z = newV3.z * R / len;
                        
                        // create the 3 new faces
                        int id0 = f.ids[0];
                        int id1 = f.ids[1];
                        int id2 = f.ids[2];
                        unsigned int idN1 = nodesGS.size();
                        for (size_t j=0; j<nodesGS.size(); ++j)
                            if (nodesGS[j] == newV1) {
                                idN1 = j;
                                break;
                            }
                        if (idN1 == nodesGS.size())
                            nodesGS.push_back(newV1);
                        
                        unsigned int idN2 = nodesGS.size();
                        for (size_t j=0; j<nodesGS.size(); ++j)
                            if (nodesGS[j] == newV2) {
                                idN2 = j;
                                break;
                                }
                        if (idN2 == nodesGS.size())
                            nodesGS.push_back(newV2);
                        
                        unsigned int idN3 = nodesGS.size();
                        for (size_t j=0; j<nodesGS.size(); ++j)
                            if (nodesGS[j] == newV3) {
                                idN3 = j;
                                break;
                            }
                        if (idN3 == nodesGS.size())
                            nodesGS.push_back(newV3);

                        Face3 f1,f2,f3,f4;
                        f1.ids[0] = id0;
                        f1.ids[1] = idN1;
                        f1.ids[2] = idN3;
                        f2.ids[0] = idN1;
                        f2.ids[1] = id1;
                        f2.ids[2] = idN2;
                        f3.ids[0] = idN2;
                        f3.ids[1] = id2;
                        f3.ids[2] = idN3;
                        f4.ids[0] = idN1;
                        f4.ids[1] = idN2;
                        f4.ids[2] = idN3;

                        facesGS.push_back(f1);
                        facesGS.push_back(f2);
                        facesGS.push_back(f3);
                        facesGS.push_back(f4);
                    }
                }

                for (size_t i=0; i<nodesGS.size(); ++i)
                    {
                        Vector3 &v = nodesGS[i];
                        double len = v.GetLength();
                        v[0] = v[0] * R / len;
                        v[1] = v[1] * R / len;
                        v[2] = v[2] * R / len;
                    }
                
                // store the result;
                //m_faces = facesGS;
                //m_nodes = nodesGS;

                // @todo : convert to our structures
                unsigned int points = nodesGS.size();
                Float3DataBlockPtr vertices =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                Float3DataBlockPtr normals =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                IDataBlockList texCoordList;            
                Float2DataBlockPtr texCoords = 
                    Float2DataBlockPtr(new DataBlock<2, float>(points));
                texCoordList.push_back(texCoords);
                Float3DataBlockPtr colors = 
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                GeometrySetPtr geom = 
                    GeometrySetPtr(new GeometrySet(vertices, 
                                                   normals, 
                                                   texCoordList,
                                                   colors));


                Vector<3, float> normal = Vector<3, float>(0, 1, 0);

                for (unsigned int i = 0; i < points; ++i){
                    Vector<3, float> v = nodesGS[i];
                    vertices->SetElement(i, v);

                    // from: http://en.wikipedia.org/wiki/UV_mapping
                    float sum = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
                    float denum = sqrt(sum);
                    float uT = v[0] / denum;
                    float vT = v[1] / denum;
                    Vector<2,float> texCoord( uT, vT );
                    texCoords->SetElement(i, texCoord);

                    Vector<3, float> n = v.GetNormalize() * inverted;
                    normals->SetElement(i, normal);

                    colors->SetElement(i, color);
                }

                unsigned int numFaces = facesGS.size();
                unsigned int numVerts = 3 * numFaces;
                unsigned int* i = new unsigned int[numVerts];
                IndicesPtr indices = IndicesPtr(new Indices(numVerts, i));
                unsigned int index = 0;
                for (unsigned int t = 0; t < numFaces; ++t){
                    for (unsigned int v = 0; v < 3; ++v){
                        i[index++] = facesGS[t].ids[v];
                    }
                }
                return MeshPtr(new Mesh(indices, TRIANGLES, 
                                        geom, MaterialPtr(new Material())));
            }

            MeshPtr CreateGeodesicSphere(float radius, unsigned int detail,
                                         bool inverted,
                                         Vector<3, float> color) {

                /*
                float size = radius;

                unsigned int d = detail + 1;
                float halfSize = size / 2;
                float unit = size / detail;
            
                unsigned int points = d * d;
                */



                MeshPtr mesh = create(detail, radius, inverted, color);
                return mesh;
            }

            MeshPtr CreateCylinder(float radius, float height,
                                             unsigned int detail,
                                             Vector<3, float> color) {
                
                unsigned int d = detail + 3;
                
                unsigned int points = 4 * d + 2;
                float radsPrSlice = 2.0 * 3.14 / float(d);

                IndicesPtr indices = IndicesPtr(new Indices(12 * d));

                Float3DataBlockPtr vertices =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                Float3DataBlockPtr normals =
                    Float3DataBlockPtr(new DataBlock<3, float>(points));
                IDataBlockList texCoords;
                Float3DataBlockPtr colors = 
                    Float3DataBlockPtr(new DataBlock<3, float>(points));

                GeometrySetPtr geom = 
                    GeometrySetPtr(new GeometrySet(vertices, 
                                                   normals, 
                                                   texCoords,
                                                   colors));
                
                unsigned int v = 0;
                unsigned int n = 0;
                unsigned int i = 0;

                // Create the lid
                vertices->SetElement(v++, Vector<3, float>(0, height / 2.0f, 0));
                vertices->SetElement(v++, Vector<3, float>(radius, height / 2.0f, 0));
                while (v < d + 1){
                    Vector<3, float> point(radius * cos((v-1) * radsPrSlice), 
                                           height / 2.0f,
                                           radius * sin((v-1) * radsPrSlice));
                    vertices->SetElement(v, point);

                    indices->GetData()[i++] = 0;
                    indices->GetData()[i++] = v - 1;
                    indices->GetData()[i++] = v;

                    v++;
                }
                indices->GetData()[i++] = 0;
                indices->GetData()[i++] = 1;
                indices->GetData()[i++] = v-1;

                for (; n < v; ++n)
                    normals->SetElement(n++, Vector<3, float>(0, 1, 0));
                
                // Create the bottom
                vertices->SetElement(v++, Vector<3, float>(0, -height / 2.0f, 0));
                vertices->SetElement(v++, Vector<3, float>(radius, -height / 2.0f, 0));
                while (v < 2 * d + 2){
                    Vector<3, float> point(radius * cos((v-d-2) * radsPrSlice), 
                                           -height / 2.0f,
                                           radius * sin((v-d-2) * radsPrSlice));
                    vertices->SetElement(v, point);

                    indices->GetData()[i++] = d + 1;
                    indices->GetData()[i++] = v - 1;
                    indices->GetData()[i++] = v;

                    v++;
                }
                indices->GetData()[i++] = d + 1;
                indices->GetData()[i++] = d + 2;
                indices->GetData()[i++] = v-1;
                
                for (; n < v; ++n)
                    normals->SetElement(n++, Vector<3, float>(0, -1, 0));

                // Create the cylinder
                Vector<3, float> topCenter(0, height / 2.0f, 0);

                vertices->SetElement(v++, Vector<3, float>(radius, height / 2.0f, 0));
                vertices->SetElement(v++, Vector<3, float>(radius, -height / 2.0f, 0));
                normals->SetElement(n++, Vector<3, float>(1, 0, 0));
                normals->SetElement(n++, Vector<3, float>(1, 0, 0));
                for (unsigned int j = 1; j < d; ++j){
                    Vector<3, float> point(radius * cos(j * radsPrSlice), 
                                           height / 2.0f,
                                           radius * sin(j * radsPrSlice));

                    Vector<3, float>  normal = (point - topCenter).GetNormalize();

                    vertices->SetElement(v++, point);
                    normals->SetElement(n++, normal);

                    indices->GetData()[i++] = v - 3;
                    indices->GetData()[i++] = v - 2;
                    indices->GetData()[i++] = v - 1;

                    point[1] = -height / 2.0f;

                    vertices->SetElement(v++, point);
                    normals->SetElement(n++, normal);

                    indices->GetData()[i++] = v - 3;
                    indices->GetData()[i++] = v - 2;
                    indices->GetData()[i++] = v - 1;
                }

                indices->GetData()[i++] = v - 2;
                indices->GetData()[i++] = v - 1;
                indices->GetData()[i++] = 2 * d + 2;
                indices->GetData()[i++] = v - 1;
                indices->GetData()[i++] = 2 * d + 2;
                indices->GetData()[i++] = 2 * d + 3;

                // Fill colors
                for (unsigned int i = 0; i < colors->GetSize(); ++i){
                    colors->SetElement(i, color);
                }

                return MeshPtr(new Mesh(indices, TRIANGLES, 
                                        geom, MaterialPtr(new Material())));
            }
                
                
        }
    }
}
