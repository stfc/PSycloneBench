# CPU CONSTANTS
# For single core of Xeon E5-1620 v2 (my desktop), as measured with 
# the Intel MKL version of linpack. This is therefore using
# 256-bit AVX instructions (SIMD)
PEAK_GFLOPS= 28.1
# An AVX instruction can do 4 d.p. FLOPs in one go so if we were
# to turn it off altogether we get a factor four reduction:
NOSIMD_PEAK_GFLOPS <- PEAK_GFLOPS/4.0

# MEM CONSTANTS
# For single core of E5-1620 v2 as measured with STREAM
# with arrays of 15M elements. Therefore, this is bandwidth to 
# main memory, not cache. Units are GB/s.
PEAK_MEM_BW= 14.5
# Using arrays of 0.25M elements I think we get bandwidth to
# L3 cache:
PEAK_L3_BW= 39.6

XMIN= 0.05
XMAX= 4.0

cpu_roof <- PEAK_GFLOPS
PEAK_BW <- max(PEAK_MEM_BW,PEAK_L3_BW)

mydata <- read.table("./roofline_haswell.dat", sep=",")
xvals <- seq(from=XMIN, to=XMAX, length.out=100)
yvals <- mydata$V2
# Width of bars to draw
barW <- 0.02

mem_roof <- function(x, peak) {
  result <- x*peak
  return(result)
}

cpu_ceiling <- function(x, peak_flops) {
    roof <- x*PEAK_BW
    vals <- rep(peak_flops, length(x)) 
    result <- pmin(roof, vals)
    return(result)
}

peak_cpu_ceiling <- function(x) cpu_ceiling(x, PEAK_GFLOPS)
nosimd_cpu_ceiling <- function(x) cpu_ceiling(x, NOSIMD_PEAK_GFLOPS)

mem_ceiling <- function(x) {
    mem_roof <- x*PEAK_MEM_BW
    peaklist <- rep(PEAK_GFLOPS, length(x))
    result <- pmin(mem_roof, peaklist)
    return(result)
}

plot(mydata$V1, mydata$V2, log="xy", xlim=c(XMIN,XMAX), ylim=c(0.01,30.0), xlab="Operational intensity (FLOPs/byte)", ylab="GFLOPS")
rect(xleft=(1.0-barW)*mydata$V1, ybottom=0.01, xright=(1.0+barW)*mydata$V1, ytop=yvals)
curve(mem_ceiling, xvals, add=TRUE, col="blue")
curve(nosimd_cpu_ceiling, xvals, add=TRUE, col="blue")
curve(peak_cpu_ceiling, xvals, add=TRUE, col="red", lty=1, lwd=2)
