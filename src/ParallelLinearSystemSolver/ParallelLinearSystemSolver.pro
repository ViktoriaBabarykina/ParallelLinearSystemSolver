TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        /usr/include/mpich

SOURCES += \
        main.c \
        matrix.c

HEADERS += \
    matrix.h
