PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /System/Library/Frameworks/Python.framework/Versions/$(PYTHON_VERSION)/include/python$(PYTHON_VERSION)

BOOST_INC = /usr/local/include
BOOST_LIB = /usr/local/Cellar/boost-python/1.60.0/lib

TARGET = ljforces

$(TARGET).so: $(TARGET).o
	clang++ -std=c++11 -shared $(TARGET).o -L$(BOOST_LIB) -lboost_python -L$(PYTHON_INCLUDE) -lpython$(PYTHON_VERSION) -o $(TARGET).so

$(TARGET).o: $(TARGET).cpp
	clang++ -std=c++11 -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cpp
