Running with 32 in horizontal and 100 in vertical for these
tests. Checking the results at the moment.

----------------

module load pgi/19.10
module load cuda/10.1

----------------

The original code running on the host -ta=host using pgi 19.10 and cuda
10.1 gives the same results if -O0, -O1 or -O are used. It gives
different results if -fast is used but -fast,-O2,-O3,-O4 and -Ofastsse
give the same results.


Performance:

-O0  81s   1.227...376E20
-O1  75s    "
-O   64s    "
-O2  64s   1.227...260E20
-fast 64s   "

So for original code -O gives as good performance as -fast but keeps
the same results as -O0.

-------------

reorder results are different. This is to be expected as the matrix
index order has been changed but the data layout has not (as we are
lazy).

However, even when the matrix order is not changed we get different
results. One reason for this is that kinner changes the order in which
the results are added up.

In the original

tmp=0.0
do i=1,ndofs
  tmp = tmp + mat(...) * v(...)
end do
result(...) = result(...) + tmp

In the ninner version

result(...) = result(...) + mat(...) * v(...)

Interestingly when the matrix vector is replaced with equivalent loops
the results change slightly for the -O[2,3,4], -fast, -fastsse results
to ...262E20 cf 260E20 in the original (but the -O[0,1] -O results
stay the same).

-------------

Just replacing the gather gives us the following results ...
1.2271898522039383E+020 (for -O0, -O1 and -O)
1.2271898522039335E+020 (-fast, -O2, -O3, -O4 and -fastsse)

However, with kinner as well it is ...
1.2271898522039380E+020 (for -O0, -O1 and -O)
1.2271898522039334E+020 (for -O2, -O3, -O4, -fast and -fastsse)

Is there any reason why kinner would change the results?  I think the
ordering of the summations changes if the thing being written to is
continuous in the vertical. Perhaps this is why? So this is another
reason for results to change.

Checking the matrix_vector_kernel_reorder1 results (without changing
the matrix index) ...

We get the same results as above 
-O0, -O1 and -O all the same as kinner above
 (1.2271898522039380E+020)
-O2, -O3, -O4, -fast, -fastsse all the same as kinner above
 (1.2271898522039334E+020)

Great. OK, we can say that the kinner code is working correctly, we
know why there are differences in results and we have baseline cpu
results to compare with on the GPU.

The baseline "correct" ninner reduction value with the matrix index
re-order is ...

1.2226489429879603E+020 (for -O0, -O1 and -O)
1.2226489429879534E+020 (for -O2, -O3, -O4, -fast and -fastsse)

-------------

Running re-order results on GPU ...

-O0, -O1 and -O same result as CPU -O
-O2, -O3, -O4, -fast and -fastsse same result as CPU -fast

Whoop!

Average of ~5 runs ignoring first as seems slower

-O0 1.68s
-O1 1.68s
-O 1.35s
-O2 1.35s
-O3 1.40s
-O4 1.40s
-fast 1.34s
-fastsse 1.34s

-------------

Checking inlined re-order version on CPU

All results the same as the not inlined version i.e. ...

1.2226489429879603E+020 (for -O0, -O1 and -O)
1.2226489429879534E+020 (for -O2, -O3, -O4, -fast and -fastsse)

Checking inlined re-order version on GPU

-O0 -O1 same
-O -O1 -O3 -O4 -fast -fastsse results vary from one run to next

Looking at why ...

If we create a temporary variable for cmap lookup ...

-O now gives same results
-O2 -O3 -O4 -fast -fastsse results vary from one run to next

Hmmm .... looks like adding a comment in the inner loop gets -O2
working most of the time!

If we create a temporary variable for writer map ...

-O2 -O3, -O4, -fast, -fastsse now give same results

Note, I've also created a temporary variable for the reader map to be
consistent, but that does not make a difference to the results.

So, it looks like a compiler bug. I should let nvidia know.

-------------

Checking performance as it looks like adding temporary vars slows the code down ...

Yes ... adding a temporary variable for the reader map slows the code down alot!!!!
So ... remove this and let NVIDIA know.

Timings with temporary cmap and writer map are ...

Average of ~5 runs ignoring first as seems slower

-O0 0.97s
-O1 1.22s
-O 0.74s
-O2 0.74s
-O3 0.76s
-O4 0.76s
-fast 0.74s
-fastsse 0.74s

For comparison 2 skylakes (not inlined) run in 0.83s.

-------------

Timings with just temporary cmap are faster but don't always give the correct reduction value ... pass to NVIDIA.

-O 0.63s - but sometimes different results.

-------------

TODO: Parallelise over reader dofs as well (use collapse)
TODO: Check naive and NVIDIA versions
TODO: Run with varying vertical
TODO: Run with varying horizontal
TODO: Check inline version on Skylake


naive gpu implementation and nvidia versions - I think they all need atomic

naive is around 51s

best nvidia is around 1.38s
All results differ and in all/most cases are not reproducible due to atomic operation.
Need to look at removing colouring.

-------------


Once I've done this I can run some plots and see how things change with different problem sizes.

I'd like to be able to test with different mappings too as some are continuous in the vertical and some are not. Of course some are discontinuous in the horizontal too and this could make a big difference.

