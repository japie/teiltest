
#ifndef _TRI_MATRIX_
#define _TRI_MATRIX_

#include <cassert>

template <class T> class TriMatrix
{

    public:
        TriMatrix ()
        {
            size = 0;
        }
        TriMatrix (int n)
        {
            size = 0;
            init(n);
        }
        void init(int n)
        {
            if (n==0) return;

            if (n==size)
            {
                for (int i=0; i<size-1; i++)
                    for (int j=0; j<i; j++)
                        matrix[i][j] = T(0.0);
            }
            else
            {
                release();
                size = n;
                matrix = new T*[size-1];
                for (int i=0; i<size-1; i++)
                    matrix[i] = new T[i+1];
            }
        }
        void release()
        {
            for (int i=0; i<size-1; i++)
                delete [](matrix[i]);
            delete []matrix;
            size = 0;
        }

        ~TriMatrix ()
        {
            release();
        }

        void write(int i, int j, T val)
        {
            assert(i<size && j<size);
            if (i<j)
            {
                int temp = i;
                i = j;
                j = temp;
            }
            matrix[i-1][j] = val;
        }

        T operator()(int i, int j)
        {
            if (i==j) return T(1.0);
            assert(i<size && j<size);
            if (i<j)
            {
                int temp = i;
                i = j;
                j = temp;
            }
            return matrix[i-1][j];
        }

    private:
        T** matrix;
        int size;

};
#endif
