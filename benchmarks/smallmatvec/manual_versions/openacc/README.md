Running with 32 in horizontal and 100 in vertical for these tests. Checking the results at the moment.

----------------

The original code running on the host -ta=host using pgi 19.9 and cuda 20 gives the same results if -O0, -O1 or -O are used. It gives different results if -fast is used but -fast and -O2 give the same results.

Performance:

-O0  81s   1.227...376E20
-O1  75s    "
-O   64s    "
-O2  64s   1.227...260E20
-fast 64s   "

So -O gives as good performance as -fast but keeps the same results as -O0.

-------------

naive gpu implementation and nvidia versions - I think they all need atomic

naive is around 51s

best nvidia is around 1.38s
All results differ and in all/most cases are not reproducible due to atomic operation.
Need to look at removing colouring.

-------------

reorder is 1.26s
reorder + inline is 0.65s
Both results differ.
Re-order gives  reproducible results.

For comparison, 2 skylakes run in 0.83s so in theory we are faster than 2 skylakes but we need to check if the results are correct.

Once I've done this I can run some plots and see how things change with different problem sizes.

I'd like to be able to test with different mappings too as some are continuous in the vertical and some are not. Of course some are discontinuous in the horizontal too and this could make a big difference.

