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

schedule=psy.invokes.get('invoke_0').schedule
#schedule.view()
new_schedule=schedule

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
