#include <cstdlib>
#include <iostream>

#include "Tsuji.h"

int main()
{
        Tsuji t;
        t.verbose();
        //t.sampling(1.67,5,60,5,0);
        BBI bb = t.getBBI(2.5,5,60,5,0);
        return EXIT_SUCCESS;



}
