# Python script which uses PSyclone to apply transformations
# before generating the PSy layer Fortran code.

def trans(psy):
    from transformations import OpenACCParallelTrans
    atrans = OpenACCParallelTrans()

    invoke = psy.invokes.get('invoke_0')
    schedule = invoke.schedule
    #schedule.view()
    new_schedule=schedule

    # Apply the OpenMP Loop transformation to *every* loop 
    # in the schedule
    for child in schedule.children:
        newschedule,memento=atrans.apply(child)
        schedule = newschedule

    invoke.schedule = newschedule
    newschedule.view()
    return psy

from parse import parse,ParseError
from psyGen import PSyFactory,GenerationError
api="gocean1.0"
filename="nemolite2d_alg.f90"
ast,invokeInfo=parse(filename,api=api,invoke_name="invoke")
psy=PSyFactory(api).create(invokeInfo)
#print psy.invokes.names

new_psy = trans(psy)

print new_psy.gen