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
of = t.get_trans_name('GOceanOpenMPLoop')

schedule=psy.invokes.get('invoke_0').schedule
#schedule.view()
new_schedule=schedule

# Apply the OpenMP transformation to *every* loop 
# in the schedule
for child in schedule.children:
    newschedule,memento=of.apply(child)
    schedule = newschedule

psy.invokes.get('invoke_0')._schedule=schedule
print psy.gen
