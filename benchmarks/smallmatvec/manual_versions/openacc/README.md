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

-O 1.26s

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

Timings with just temporary cmap are faster but don't always give the
correct reduction value ... pass to NVIDIA.

-O 0.63s - but sometimes different results.

-------------

Parallelising over reader dofs as well using collapse(2) just makes
the code run more slowly (~1.45s) even when the matrix 2nd and 3rd
dimensions are swapped. I presume this is to do with the sparse read
to the rhs vector.

I guess there might be a regime where we need more parallelism so this
would be beneficial?

-------------

Ran reorder and inlined reorder with varying vertical
Ran reorder and inlined with varying horizontal
Results are in results/*
Comparison results with CPU are in ../results/*

-------------

Manually made matrix the original (bad) way round i.e. we now have
kinner and kinner inlined without re-order.
Results are in results/*

-------------

Looking at naive and nvidia versions

naive version needs atomic as we are writing to a field that is
continuous in the vertical. If it were not we would not need
atomic. We also need to privatise the two arrays gather arrays.

time                checksum
64.76465797424316s  1.2213850309400101E+020

nvidia1
3.580666065216064   1.3805756942507162E+020

nvidia2
3.909955978393555   1.3805756942507160E+020

nvidia3 - core dump

nvidia4 - core dump

-------------

Running naive and nvidia versions serially on cpu with -O0 to check for
correct coding and fixing any bugs ...

naive version: 1.2271898522039376E+020 same as original code
nvidia1 version: 1.2271898522039376E+020 same as original code
nvidia2 version: 1.2271898522039376E+020 same as original code
nvidia3 version: 1.2271898522039376E+020 same as original code
nvidia4 version: 1.2271898522039376E+020 same as original code

Whoop!

Now running on GPU ...

naive version:
-O0 64.90566587448120 1.2213850309400101E+020
-O1 64.89450097084045 1.2213850309400101E+020
-O 64.77101993560791 1.2213850309400101E+020
-O2 64.77753400802612 1.2213850309400125E+020
-O3 64.93016695976257 1.2213850309400125E+020
-O4 64.91746497154236 1.2213850309400125E+020
-fast 64.76277613639832 1.2213850309400125E+020
-fastsse 64.78245496749878 1.2213850309400125E+020

*** Naive checksum looks dodgy but can't see problem unless it is private arrays?

nvidia1
-O0 3.501081943511963 1.2271898522039380E+020
-O1 3.455302953720093 1.2271898522039380E+020
-O 3.269083023071289 1.2271898522039380E+020
-O2 3.306182861328125 1.2271898522039262E+020
-O3 3.783846855163574 1.2271898522039262E+020
-O4 3.792832136154175 1.2271898522039262E+020
-fast 3.264442920684814 1.2271898522039262E+020
-fastsse 3.256054878234863 1.2271898522039262E+020

nvidia2
-O0 4.415854930877686 1.2271898522039376E+020
-O1 4.449892044067383 1.2271898522039376E+020
-O 3.573765993118286 1.2271898522039376E+020
-O2 3.570590972900391 1.2271898522039262E+020
-O3 4.251459121704102 1.2271898522039262E+020
-O4 4.239311933517456 1.2271898522039262E+020
-fast 3.563273906707764 1.2271898522039262E+020
-fastsse 3.559406995773315 1.2271898522039262E+020

nvidia3
-O0 1.735004901885986 1.2271898522039376E+020
-O1 1.714452028274536 1.2271898522039376E+020
-O 1.440443992614746 1.2271898522039376E+020
-O2 1.448733091354370 1.2271898522039262E+020
-O3 1.470628976821899 1.2271898522039262E+020
-O4 1.459609985351563 1.2271898522039262E+020
-fast 1.448508024215698 1.2271898522039262E+020
-fastsse 1.456532955169678 1.2271898522039262E+020

nvidia4
-O0 1.082563161849976 1.2271898522039376E+020
-O1 1.081111907958984 1.2271898522039380E+020
-O 0.9792480468750000 1.2271898522039380E+020
-O2 0.9680049419403076 1.2271898522039262E+020
-O3 0.9599709510803223 1.2271898522039262E+020
-O4 0.9638509750366211 1.2271898522039262E+020
-fast 0.9689779281616211 1.2271898522039262E+020
-fastsse 0.9708011150360107 1.2271898522039262E+020


So naive results look dodgy but are slow so we don't really care.

nvidia4 gives good performance and is the fastest
So let's run with this version and get results to plot ...
Using -O

-------------

Created an atomic GPU version with matrix reorder and collapse colour loops

reorder3 results are:
-O0 0.6700620651245117 1.2226489429879602E+020
-O1 0.6696147918701172 1.2226489429879602E+020
-O 0.6656098365783691 1.2226489429879603E+020
-O2 0.6654019355773926 1.2226489429879536E+020
-O3 0.6650259494781494 1.2226489429879536E+020
-O4 0.6659181118011475 1.2226489429879536E+020
-fast 0.6659500598907471 1.2226489429879536E+020
-fastsse 0.6659710407257080 1.2226489429879536E+020

Notice the results are always the same. This problem must be
discontinuous in the horizontal.

Running with this version to get results to plot ...
Using -O

-------------

UP TO HERE

After this I will try to get this onto master and leave other changes for future PRs

TODO: is colouring really needed in this example? - seems not.
TODO: inline version on Skylake
TODO: restructure setup code
TODO: MPI + OpenMP on Skylake
TODO: remove colouring and run GPU results again
TODO: revisit locking with skylake

-------------

I'd like to be able to test with different mappings as some are
continuous in the vertical and some are not. Also, some are
discontinuous in the horizontal too and this could make a difference
(no locking or colouring required).

