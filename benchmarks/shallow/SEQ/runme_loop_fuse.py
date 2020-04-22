from parse import parse,ParseError
from psyGen import PSyFactory,GenerationError
#from algGen import Alg
api="gocean"
filename="shallow_gocean.f90"
ast,invokeInfo=parse(filename,api=api,invoke_name="invoke")
psy=PSyFactory(api).create(invokeInfo)
print psy.gen
#alg=Alg(ast,psy)

print psy.invokes.names
schedule=psy.invokes.get('invoke_0').schedule
schedule.view()

from psyGen import TransInfo
t=TransInfo()
print t.list
#lf=t.get_trans_name('DoubleLoopFuse')
lf=t.get_trans_name('LoopFuse')

newschedule,memento=lf.apply(schedule.children[0],schedule.children[1])
#newschedule,memento=lf.apply(schedule.children[0].children[0].children[0],schedule.children[1].children[0].children[0])
newschedule.view()
#psy.invokes.get('invoke_0')._schedule=newschedule
#print psy.gen
