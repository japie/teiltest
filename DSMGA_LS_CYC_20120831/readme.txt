/*******************************/
/Author:        ChungHsiang     /
/Last modified: 2010.10.16      /
/*******************************/

(1) Replace all files which own the same name as those are contained in "modified.zip" 

(2) type "make clean" to remove currently exist dependence

(3) type "make dataSampling" to generate executable "dataSampling"

(4) type "make" to generate executable DSMGA

/*********************dataSampling********************/
./dataSampling ell numGroup std
ell:      length of the chromosome;
numGroup: ell/interval of the mean;recommended value = ell / 3 in Monetomo set
std =     standard deviation = sigma;

The data will be generated as "NormalSet/NSet_ell_numGroup_std.dat"

/*************************DSMGA**********************/
./DSMGA ell nInitial SP pc pm maxGen maxFe repeat display rand_seed std
-> an new argument "std" is added to the commmand

#interval is defined in global.h
#if (std != -1)
 ->load the test data specified by ell and std(if exist)
 ->getBBIMone(data) is called after object DSMGA is constructed

