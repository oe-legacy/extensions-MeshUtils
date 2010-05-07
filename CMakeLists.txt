SET ( EXTENSION_NAME "Extensions_MeshUtils")

# Create the extension library
ADD_LIBRARY( ${EXTENSION_NAME}
  Utils/MeshCreator.h
  Utils/MeshCreator.cpp
  Utils/MeshTransformer.h
  Utils/MeshTransformer.cpp
  Utils/MeshUtils.h
  Utils/MeshUtils.cpp
)

TARGET_LINK_LIBRARIES( ${EXTENSION_NAME}
  OpenEngine_Resources
  OpenEngine_Geometry
)
