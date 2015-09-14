# Python script which uses PSyclone to apply transformations
# before generating the PSy layer Fortran code.
from parse import parse,ParseError
from psyGen import PSyFactory,GenerationError
api="gocean1.0"
filename="nemolite2d_alg.f90"
ast,invokeInfo=parse(filename,api=api,invoke_name="invoke")
psy=PSyFactory(api).create(invokeInfo)
#print psy.invokes.names

from psyGen import TransInfo
t=TransInfo()
ltrans = t.get_trans_name('GOceanOMPLoopTrans')
rtrans = t.get_trans_name('OMPParallelTrans')
ftrans = t.get_trans_name('GOceanLoopFuse')

schedule=psy.invokes.get('invoke_0').schedule
#schedule.view()
new_schedule=schedule

# Fuse the three field-copy operations
# First the outer loops...
new_schedule, memento = ftrans.apply(schedule.children[8],
                                     schedule.children[9])
new_schedule, memento = ftrans.apply(new_schedule.children[8],
                                     new_schedule.children[9])
# Then the inner loops...
new_schedule, memento = ftrans.apply(schedule.children[8].children[0],
                                     schedule.children[8].children[1])
new_schedule, memento = ftrans.apply(schedule.children[8].children[0],
                                     schedule.children[8].children[1])

# Apply the OpenMP Loop transformation to *every* loop 
# in the schedule
for child in schedule.children:
    newschedule,memento=ltrans.apply(child)
    schedule = newschedule


# Enclose all of these loops within a single OpenMP
# PARALLEL region
newschedule,memento = rtrans.apply(schedule.children)

psy.invokes.get('invoke_0')._schedule=schedule
print psy.gen
