#  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick
#  This code is licensed under the MIT License.

# https://stackoverflow.com/questions/40621451/makefile-automatically-compile-all-c-files-keeping-o-files-in-separate-folde

SRC=src
LIB=src
BIN=bin
OBJ=bin
SO=bin
TEST=test
EXP=experiments
PYDIR=python

DEBUG=
# https://stackoverflow.com/questions/3375697/what-are-the-useful-gcc-flags-for-c
SCIP=
PYTHON=t

ifdef DEBUG
CPPFLAGS=-g3 -fPIC -Wall -O0 #-Wextra -Wshadow
LDFLAGS=-g3
else
CPPFLAGS=-fPIC -Wall -O2 -I/opt/homebrew/include
LDFLAGS=-L/opt/homebrew/lib
endif

ifdef SCIP
LDLIBS=-lflint -lboost_iostreams -lboost_system -lboost_filesystem -lScipPP -lscip
else
LDLIBS=-lflint -lboost_iostreams
endif

ifdef PYTHON
PYLIBS=$(shell python3 -m pybind11 --includes)
PYEXT=$(shell python3-config --extension-suffix)
endif

CC=g++

#CPPFLAGS=-g -Wall -Werror -Wextra -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wstrict-overflow=5 -Wwrite-strings -Waggregate-return -Wcast-qual -Wswitch-default -Wswitch-enum -Wconversion -Wunreachable-code -Wformat=2 -O2


LIB_SRC=$(wildcard $(SRC)/*.cpp)
LIB_OBJ=$(patsubst $(SRC)/%.cpp, $(OBJ)/lib_%.o, $(LIB_SRC))
LIB_HPP=$(wildcard $(LIB)/*.hpp)

TEST_SRC=$(wildcard $(TEST)/*.cpp)
TEST_OBJ=$(patsubst $(TEST)/%.cpp, $(OBJ)/test_%.o, $(TEST_SRC))
TEST_BIN=$(patsubst $(TEST)/%.cpp, $(BIN)/test_%.bin, $(TEST_SRC))

EXP_SRC=$(wildcard $(EXP)/*.cpp)
EXP_HPP=$(wildcard $(EXP)/*.hpp)
EXP_OBJ=$(patsubst $(EXP)/%.cpp, $(OBJ)/exp_%.o, $(EXP_SRC))
EXP_BIN=$(patsubst $(EXP)/%.cpp, $(BIN)/exp_%.bin, $(EXP_SRC))

ifdef PYTHON
PY_SRC=$(wildcard $(PYDIR)/*.cpp)
PY_HPP=$(wildcard $(PYDIR)/*.hpp)
PY_OBJ=$(patsubst $(PYDIR)/%.cpp, $(OBJ)/py_%.o, $(PY_SRC))
PY_SO=$(patsubst $(PYDIR)/%.cpp, $(SO)/%$(PYEXT), $(PY_SRC))
endif

ALL_OBJ=$(LIB_OBJ) $(TEST_OBJ) $(EXP_OBJ)
ALL_BIN=$(TEST_BIN) $(EXP_BIN) $(PY_SO)

all: $(ALL_BIN)

$(BIN)/test_%.bin: $(OBJ)/test_%.o $(LIB_OBJ)
	$(CC) $(LDFLAGS) -o $@ $< $(LIB_OBJ) $(LDLIBS)

$(BIN)/exp_%.bin: $(OBJ)/exp_%.o $(LIB_OBJ)
	$(CC) $(LDFLAGS) -o $@ $< $(LIB_OBJ) $(LDLIBS)

$(OBJ)/lib_%.o: $(SRC)/%.cpp $(LIB_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@ 

$(OBJ)/test_%.o: $(TEST)/%.cpp $(LIB_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@ 

$(OBJ)/exp_%.o: $(EXP)/%.cpp $(LIB_HPP) $(EXP_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) -c $< -o $@

$(OBJ)/py_%.o: $(PYDIR)/%.cpp $(LIB_HPP) $(PY_HPP)
	$(CC) $(CPPFLAGS) -I$(LIB) $(PYLIBS) -c $< -o $@

$(OBJ)/%$(PYEXT): $(OBJ)/py_%.o $(LIB_OBJ)
	$(CC) $(LDFLAGS) -shared -o $@ $< $(LIB_OBJ) $(LDLIBS)


.PHONY clean:
	rm -f $(SO)/*$(PYEXT)
	rm -f $(OBJ)/*.o
	rm -f $(BIN)/*.bin
