/*
 * Simple parallel quicksort, check wikipedia
 * Shamelessly hacked from glibc, thereby GPL2
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience;
 * Vol. 23 (11), 1249-1265, 1993.
 */

#include "globals.h"
#include <gsl/gsl_heapsort.h>

#define PARALLEL_THRES_QSORT 15000 // use serial sort below this limit
#define PARALLEL_THRES_HEAPSORT 15000

#define INSERT_THRES 8 // insertion sort threshold

#define SWAP(a, b, size)             \
    do {                            \
        size_t __size = (size);        \
        char *__a = (a), *__b = (b);    \
        do {                            \
            char __tmp = *__a;            \
            *__a++ = *__b;                \
            *__b++ = __tmp;                \
        } while (--__size > 0);        \
    } while (0)

#define SWAP_SIZE_T(a,b)                \
    do {                                \
        size_t *__a = (a), *__b = (b);    \
        size_t __tmp = *__a;            \
        *__a = *__b;                    \
        *__b = __tmp;                    \
    } while (0)


#define COMPARE_DATA(a,b,size) ((*cmp) ((void *)(data + *a * size), \
            (void *)(data + *b * size))) // compare with non-permutated data

#define STACK_SIZE (CHAR_BIT * sizeof(size_t))

/*
 * A shared stack keeps the partition pointers. We assign one thread per
 * partition until there is a partition for every thread.
 * Then we continue serial on every partition. There is little
 * need for synchronisation, which is done intrinsically by the for(j) loop.
 * Once the array is presorted the standard qsort lib routine excels
 */

static struct SharedStackDataChar { // work queue, holding partitions
    char *lo;
    char *hi;
} shared_stack_char[STACK_SIZE] = {{NULL,NULL}};

void Qsort(const int nThreads, void *const data_ptr, int nData, size_t size,
        int (*cmp) (const void *, const void *))
{
    if (nData < PARALLEL_THRES_QSORT || nThreads == 1) {

#pragma omp single
        qsort(data_ptr, nData, size, cmp);

        return ;
    }

    /* partition number and size */
    const int desNumPar = max(1, 2*min(nThreads/2, nData/INSERT_THRES/2));

#pragma omp single  // initial stack node is just the whole array
    {

        shared_stack_char[0].lo = (char *)data_ptr;
        shared_stack_char[0].hi = &((char *)data_ptr)[size * (nData - 1)];

    } // omp single

    int delta = desNumPar << 1;

    while (delta > 2) {

        delta >>= 1;

#pragma omp for schedule (static,1)
        for (int j = 0; j < desNumPar; j += delta) {

            char *lo = shared_stack_char[j].lo; // pop from stack
            char *hi = shared_stack_char[j].hi;

            /* Find pivot element from median and sort the three.
             * That helps to prevent the n^2 worst case */
            char *mid = lo + size * ((hi - lo) / size >> 1); // pivot

            if ( (*cmp) ((void *) mid, (void *) lo) < 0)
                SWAP(mid, lo, size);

            if ( (*cmp) ((void *) hi, (void *) mid) < 0)
                SWAP (mid, hi, size);
            else
                goto jump_over;

            if ( (*cmp) ((void *) mid, (void *) lo) < 0)
                SWAP (mid, lo, size);

jump_over:;

          char *left_ptr  = lo + size;
          char *right_ptr = hi - size;

          /* now put all larger/smaller than the pivot on the right/left */
          do {

              while ((*cmp)((void *)left_ptr, (void *)mid) < 0)
                  left_ptr += size;

              while ((*cmp)((void *)mid,(void *)right_ptr) < 0)
                  right_ptr -= size;

              if (left_ptr < right_ptr) {

                  SWAP (left_ptr, right_ptr, size);

                  if (mid == left_ptr)
                      mid = right_ptr;
                  else if (mid == right_ptr)
                      mid = left_ptr;

                  left_ptr += size;
                  right_ptr -= size;

              } else if (left_ptr == right_ptr) {

                  left_ptr += size;
                  right_ptr -= size;

                  break;
              }

          } while (left_ptr <= right_ptr);

          /* Push next iterations / partitions to the shared stack */

          shared_stack_char[j].lo = lo;
          shared_stack_char[j].hi = right_ptr;

          shared_stack_char[j + (delta>>1)].lo = left_ptr;
          shared_stack_char[j + (delta>>1)].hi = hi;
        } // j

#pragma omp barrier

    } // while

#pragma omp for schedule(static,1)
    for (int i = 0; i < desNumPar; i++) { // serial sort on subpartitions

        const char *hi = shared_stack_char[i].hi;
        const char *lo = shared_stack_char[i].lo;

        size_t partition_size = (hi - lo + 1) / size;

        qsort((void *) lo, partition_size, size, cmp); // hard to beat
    } // i

    return ;
}

/*
 * This is an OpenMP parallel external sort.
 * We use the same algorithm as above to divide and conquer.
 * In the first stage we sort 1 partition per thread.
 * The work is inserted in the shared stack with a delta that shrinks
 * from N partitions / 2 to 1. Some partitions can be empty, when a
 * previous partition is too small to split.
 * The second stage is a threaded qsort on the subpartitions down
 * to INSERT_THRES followed by insertion sort on the nearly
 * ordered subpartition.
 * Here we are swapping ONLY the permutation array *perm and
 * are comparing ONLY to the data array indexed by *perm
 */

static struct SharedStackDataSizeT {
    size_t *lo;
    size_t *hi;
} shared_stack_sizet[STACK_SIZE] = {{0,0}};

void Qsort_Index(const int nThreads, size_t *perm, void *const data,
        const int nData, const size_t datasize,
        int (*cmp) (const void *, const void *))
{
    //if (nData < PARALLEL_THRES_HEAPSORT || nThreads == 1) {

#pragma omp single
    gsl_heapsort_index(perm, data, nData, datasize, cmp);

    return ;
    //}

    const int desNumPar = min(nThreads, floor(nData/INSERT_THRES));

#pragma omp for
    for (size_t i = 0; i < nData; i++ )
        perm[i] = i;

#pragma omp single
    {

        shared_stack_sizet[0].lo = perm;
        shared_stack_sizet[0].hi = &perm[nData - 1];

    } // omp single

    int delta = desNumPar << 1; //twice the distance of entries in shared stack

    /* First stage: subpartitions, roughly NThreads */

    while (delta > 2) {

        delta >>= 1; // we compensated for this before

#pragma omp flush

#pragma omp for schedule(static,1)
        for (int i = 0; i < desNumPar; i += delta) {

            size_t *lo = shared_stack_sizet[i].lo;
            size_t *hi = shared_stack_sizet[i].hi;

            if (hi - lo < 2)
                continue;

            size_t *mid = lo + ((hi - lo) >> 1); // pivot & presort

            if (COMPARE_DATA(mid,lo,datasize) < 0)
                SWAP_SIZE_T (mid, lo);

            if (COMPARE_DATA(hi,mid,datasize) < 0)
                SWAP_SIZE_T (hi, mid);
            else
                goto jump_over;

            if (COMPARE_DATA(mid,lo,datasize) < 0)
                SWAP_SIZE_T (mid, lo);

jump_over:;

          size_t *left  = lo + 1;
          size_t *right = hi - 1;

          do { // collapse the walls

              while (COMPARE_DATA(left,mid,datasize) < 0)
                  left++;

              while (COMPARE_DATA(mid,right,datasize) < 0)
                  right--;

              if (left < right) {

                  SWAP_SIZE_T (left, right);

                  if (mid == left)
                      mid = right;
                  else if (mid == right)
                      mid = left;

                  left++;
                  right--;

              } else if (left == right) {

                  left++;

                  break;
              }

          } while (left <= right);

          shared_stack_sizet[i].lo = lo; // Push to replace old position
          shared_stack_sizet[i].hi = right;

          shared_stack_sizet[i + (delta>>1)].lo = left; //Push for next
          shared_stack_sizet[i + (delta>>1)].hi = hi;
        } // j

    } // while

    /* Second stage, every thread: Qsort on subpartition *beg to *end */

#pragma omp flush

    size_t *beg = shared_stack_sizet[Omp.ThreadID].lo;
    size_t *end = shared_stack_sizet[Omp.ThreadID].hi;

    size_t partition_size = end - beg + 1;

    if (partition_size < INSERT_THRES)
        goto insertion_sort;

    struct SharedStackDataSizeT stack[STACK_SIZE] = { {beg,end} };

    int next = 0;

    for (;;)  {

        size_t *lo = stack[next].lo;  // pop from private stack
        size_t *hi = stack[next--].hi;

        size_t *mid = lo + ((hi - lo) >> 1);

        if (COMPARE_DATA(mid,lo,datasize) < 0)
            SWAP_SIZE_T (mid, lo);

        if (COMPARE_DATA(hi,mid,datasize) < 0)
            SWAP_SIZE_T (mid, hi);
        else
            goto hop_over;

        if (COMPARE_DATA(mid,lo,datasize) < 0)
            SWAP_SIZE_T (mid, lo);

hop_over:;

         size_t *left  = lo + 1;
         size_t *right = hi - 1;

         do { // collapse the walls

             while (COMPARE_DATA(left,mid,datasize) < 0)
                 left++;

             while (COMPARE_DATA(mid,right,datasize) < 0)
                 right--;

             if (left < right) {

                 SWAP_SIZE_T (left, right);

                 if (mid == left)
                     mid = right;
                 else if (mid == right)
                     mid = left;

                 left++;
                 right--;

             } else if (left == right) {

                 left++;
                 right--;

                 break;
             }

         } while (left <= right);

         if (right - lo < hi - left) { // push smaller partition first

             if (right - lo + 1 > INSERT_THRES) { // push

                 stack[++next].lo = lo;
                 stack[next].hi = right;
             }

             if (hi - left + 1 > INSERT_THRES) {

                 stack[++next].lo = left;
                 stack[next].hi = hi;
             }

         } else {

             if (hi - left + 1 > INSERT_THRES) {

                 stack[++next].lo = left;
                 stack[next].hi = hi;
             }

             if (right - lo + 1 > INSERT_THRES) { // push

                 stack[++next].lo = lo;
                 stack[next].hi = right;
             }
         }

         if (next < 0) // stack empty
             break;

    } // for()

    /* Third stage, every thread: insertion sort on its partition */

insertion_sort:;

               size_t *trail = beg;

               const size_t *runMax = min(end, beg + INSERT_THRES);

               for (size_t *run = beg+1; run <= runMax; run++) // smallest element first
                   if (COMPARE_DATA(run, trail, datasize) < 0)
                       trail = run;

               if (trail != beg)
                   SWAP_SIZE_T(trail, beg);

               size_t *run = beg + 1;

               while (++run <= end) { // insertion sort left to right

                   trail = run - 1;

                   while (COMPARE_DATA(run, trail, datasize) < 0) // find correct place
                       trail--;

                   trail++;

                   if (trail != run)  { // now move block one forward

                       size_t save = *run;

                       size_t *hi = run;
                       size_t *lo = hi - 1;

                       while (lo >= trail)
                           *hi-- = *lo--;

                       *hi = save; // sort in *run element
                   }
               } // while

#pragma omp barrier

               return;
}
