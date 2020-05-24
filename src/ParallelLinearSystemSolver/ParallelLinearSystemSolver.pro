TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        /usr/include/mpich

SOURCES += \
        common_math.c \
        conjugate_gradient.c \
        conjugate_gradient_parallel.c \
        debugging.c \
        main.c \
        matrix.c

HEADERS += \
    common_math.h \
    conjugate_gradient.h \
    conjugate_gradient_parallel.h \
    debugging.h \
    matrix.h
