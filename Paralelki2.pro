TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/lib/gcc/x86_64-linux-gnu/7/include

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
LIBS += -lgomp

SOURCES += \
        BasicOperation.cpp \
        CommandLineArg.cpp \
        Grid.cpp \
        MPI_Exchange.cpp \
        MatrixCSR.cpp \
        Solver.cpp \
        Vector.cpp \
        main.cpp

HEADERS += \
    BasicOperation.h \
    CommandLineArg.h \
    Grid.h \
    MPI_Exchange.h \
    MatrixCSR.h \
    Solver.h \
    Vector.h
