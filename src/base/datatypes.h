#ifndef DATATYPES_H
#define DATATYPES_H


template<class T>
struct DataType
{
    T data;

    DataType() {}
    DataType(T _data) : data(_data) {}
};

#endif // DATATYPES_H
