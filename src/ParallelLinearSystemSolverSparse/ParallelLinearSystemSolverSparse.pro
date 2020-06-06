TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        /usr/lib/gcc/x86_64-linux-gnu/7/include \
        /usr/include/mpich

SOURCES += \
    SparseMatrix.cpp \
    Vector.cpp \
        conj_grad_solve.cpp \
        IO.cpp \
        main.cpp \

HEADERS += \
    Matrix.hpp \
    SparseMatrix.hpp \
    Vector.hpp \
    	conj_grad_solve.hpp \
        IO.hpp \
