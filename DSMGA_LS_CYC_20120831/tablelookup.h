
#ifndef _TABLE_LOOK_UP_H
#define _TABLE_LOOK_UP_H

typedef enum TableLookUpType {
    hIFF,
    hXOR,
    hTrap
};

class TableLookUp {

public:
    TableLookUp();
    ~TableLookUp();

    void readTable(TableLookUpType type);
    void createTable(TableLookUpType type);
    
    bool lookUpBool(int);
    int lookUpInt(int);
    double lookUpDouble(int);

private:
    
    bool *tableBool;
    int *tableInt;
    double *tableDouble;
};

#endif

