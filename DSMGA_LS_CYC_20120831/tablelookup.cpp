
#include <stdio.h>
#include "chromosome.h"
#include "tablelookup.h"

TableLookUp::TableLookUp()
{
    tableInt = new int [(1<<16)];
    tableBool = new bool [(1<<16)];
    tableDouble = new double [(1<<9)];
}
    
TableLookUp::~TableLookUp()
{
    delete []tableInt;
    delete []tableBool;
    delete []tableDouble;
}

void TableLookUp::createTable(TableLookUpType type)
{
    int i;
    FILE *fp;
    
    switch (type) {
	case hIFF:
	    for (i=0; i<(1<<16); i++) {
		Chromosome ch(16);
		ch.gene[0] = i;
		tableInt[i] = ch.hIFFBigH(0, 16, false);
		tableBool[i] = ch.hIFFSmallH(0, 16, false);
	    }
	    fp = fopen("hiff-bigh.table", "wb");
	    fwrite(tableInt, sizeof(int), (1<<16), fp);
	    fclose(fp);
	    fp = fopen("hiff-smallh.table", "wb");
	    fwrite(tableBool, sizeof(bool), (1<<16), fp);
	    fclose(fp);
	    break;
	    
	case hXOR:
            for (i=0; i<(1<<16); i++) {
                Chromosome ch(16);
                ch.gene[0] = i;
                tableInt[i] = ch.hXORBigH(0, 16, false);
                tableBool[i] = ch.hXORSmallH(0, 16, false);
            }
            fp = fopen("hxor-bigh.table", "wb");
            fwrite(tableInt, sizeof(int), (1<<16), fp);
            fclose(fp);
            fp = fopen("hxor-smallh.table", "wb");
            fwrite(tableBool, sizeof(bool), (1<<16), fp);
            fclose(fp);
	    break;

	case hTrap:
	    break;
    }
			
}

void TableLookUp::readTable(TableLookUpType type)
{
    FILE *fp1;
    FILE *fp2;

    
    switch (type) {
	case hIFF:
	    fp1 = fopen("hiff-smallh.table", "rb");
	    fp2 = fopen("hiff-bigh.table", "rb");
	    if (fp1 == NULL || fp2 == NULL) {
		if (fp1 != NULL) fclose(fp1);
		if (fp2 != NULL) fclose(fp2);
		createTable(type);
	    }
	    else {
		fread(tableBool, sizeof(bool), (1<<16), fp1);
		fread(tableInt, sizeof(int), (1<<16), fp2);
		fclose(fp1);
		fclose(fp2);
	    }
		
	    break;

	case hXOR:
	    fp1 = fopen("hxor-smallh.table", "rb");
	    fp2 = fopen("hxor-bigh.table", "rb");
	    if (fp1 == NULL || fp2 == NULL) {
		if (fp1 != NULL) fclose(fp1);
		if (fp2 != NULL) fclose(fp2);
		createTable(type);
	    }
	    else {
		fread(tableBool, sizeof(bool), (1<<16), fp1);
		fread(tableInt, sizeof(int), (1<<16), fp2);
		fclose(fp1);
		fclose(fp2);
	    }
	    break;

	case hTrap:
	    break;
    }
			
}

bool TableLookUp::lookUpBool(int index)
{
    return tableBool[index];
}

int TableLookUp::lookUpInt(int index)
{
    return tableInt[index];
}

double TableLookUp::lookUpDouble(int index)
{
    return tableDouble[index];
}
    
