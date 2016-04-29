TEMPLATE = app
CONFIG += qt warn_on release embed_manifest_exe 
CONFIG -= app_bundle  # don't build Mac app bundle, for now
QT += opengl xml widgets
FORMS += edge.ui param.ui measure.ui quickselect.ui 
SOURCES += edgegui.cpp 3dview.cpp marchingcubes.cpp geometry.cpp mesh.cpp volume.cpp edgeconfig.cpp
SOURCES += edge.cpp tribox3.cpp procthreads.cpp
HEADERS += edgegui.hpp 3dview.hpp marchingcubes.hpp geometry.hpp mesh.hpp volume.hpp edgeconfig.hpp
HEADERS += edge.hpp procthreads.hpp
HEADERS += util.hpp
QMAKE_CXXFLAGS += -msse2 -ftree-vectorize 
QMAKE_CXXFLAGS += -I/usr/local/include 
#QMAKE_CXXFLAGS += -msse2 -ftree-vectorize -ftree-vectorizer-verbose=3
#macx {
# Note that for Linux it's best to compile and then install tiff-3.9.5 into /usr/local
# Make sure you use ./configure --disable-jpeg when building libtiff.
## Use statist libtiff when compiling binary for mac users.
#LIBS += -L/usr/local/lib /usr/local/lib/libjpeg.a /usr/local/lib/libtiff.a  -lz
QMAKE_CXXFLAGS += -I/usr/local/include 
LIBS += -L/usr/local/lib -ltiff  -lz 
#}
#unix {
#LIBS += -ltiff
#}
win32 {
   QMAKE_CFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
   QMAKE_CXXFLAGS += "/D_HAS_ITERATOR_DEBUGGING=0 /D_SECURE_SCL=0"
   INCLUDEPATH += "C:/tiff-4.0.3/libtiff"
   LIBS += "C:/tiff-4.0.3/libtiff/libtiff.lib"
}
