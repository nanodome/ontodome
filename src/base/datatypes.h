#ifndef DATATYPES_H
#define DATATYPES_H

/// Generic data type base template class.
/// It provides the data entry to be filled with any data
/// such as double, float, string etc.

template<class T>
struct DataType
{
    T data;

    DataType() {}
    DataType(T _data) : data(_data) {}
};

#endif // DATATYPES_H
