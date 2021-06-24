import numpy as np
import dimod
from dwave.system import DWaveSampler, AutoEmbeddingComposite
from dwave.inspector import show
from dwave.cloud import Client
import random

#Create problem instance

import numpy as np
import timeit

def planted_solution(network,loop_min=4, loop_max=100, number_of_loops=10,scale=1):
    
    n=len(network)
    solution_str=np.random.randint(2**n)
    solution_str=format(solution_str, '0'+str(n)+'b')
    solution=list(solution_str)
    loop_number=0

    ising=np.zeros((n,n))

    while loop_number<number_of_loops:
        path=[]
        ising_local=np.zeros((n,n))
        #Choose initial value
        path.append(np.random.choice(list(range(len(network)))))

        #Create path until it collides with it self
        while len(path)==len(set(path)):
            path.append(np.random.choice(network[str(path[-1])]))

        #Remove the tails and allow for length selction
        start=[k for k, y in enumerate(path) if y==path[-1]]
        loop=[x for ii, x in enumerate(path) if ii>= start[0]]

        loop_record=[]

        if loop_min<len(loop)<loop_max and not(loop in loop_record):
            #Set Ising parameters
            loop_record=loop_record.append(loop)
            for k, connection in enumerate(loop[:-1]):
                con_bit=loop[k+1]
                if solution[connection]==solution[con_bit]:
                    ising_local[connection][con_bit]=-1
                else:
                    ising_local[connection][con_bit]=1
            #Create frustration
            flip=np.random.randint(len(loop)-1)
            ising_local[loop[flip]][loop[flip+1]]=ising_local[loop[flip]][loop[flip+1]]*-1

            #Update the ising model and counter
            ising=ising+ising_local
            loop_number=loop_number+1

    ising=ising+np.transpose(ising)

    if scale:
        maxValue=np.amax(np.absolute(ising))
        ising=ising/maxValue*scale

    


    return{'solution':solution_str,'ising_model':ising}

#Function to convert matrix representation to a dictionary:
def adjgraph(adj):
    adj=np.array(adj)
    uptri=np.triu(adj,1)
    
    h={int(qubit):singlefield for qubit,singlefield in enumerate(np.diagonal(adj))}
    j={(qubit1,qubit2):interaction for qubit1,conn in enumerate(uptri) for qubit2,interaction in enumerate(conn) if interaction!=0}

    return h,j

#Function to calculate energy of a string
def energy(string, adj):
    n=len(adj)
    h=np.diag(adj)
    adj=np.array(adj/2)-np.diag(h)/2

    bits=list(string)
    zeig=np.array([1-2*int(x) for x in bits])
    value=np.dot(zeig,np.dot(adj,zeig))+np.dot(np.array(h),zeig)
    

    return {'energy':value,'string':string}




#Connections in the unit cell
unit_cell={'0':[4,5,6,7],'1':[4,5,6,7],'2':[4,5,6,7],'3':[4,5,6,7],'4':[0,1,2,3],'5':[0,1,2,3],'6':[0,1,2,3],'7':[0,1,2,3]}


chimera16={'0':[4,5,6,7],'1':[4,5,6,7],'2':[4,5,6,7],'3':[4,5,6,7],'4':[0,1,2,3,12],'5':[0,1,2,3,13],'6':[0,1,2,3,14],'7':[0,1,2,3,15],
'8':[12,13,14,15],'9':[12,13,14,15],'10':[12,13,14,15],'11':[12,13,14,15],'12':[8,9,10,11,4],'13':[8,9,10,11,5],'14':[8,9,10,11,6],'15':[8,9,10,11,7]}

chimera24={'0':[4,5,6,7,16],'1':[4,5,6,7,17],'2':[4,5,6,7,18],'3':[4,5,6,7,19],'4':[0,1,2,3,12],'5':[0,1,2,3,13],'6':[0,1,2,3,14],'7':[0,1,2,3,15],
'8':[12,13,14,15],'9':[12,13,14,15],'10':[12,13,14,15],'11':[12,13,14,15],'12':[8,9,10,11,4],'13':[8,9,10,11,5],'14':[8,9,10,11,6],'15':[8,9,10,11,7],
'16':[20,21,22,23,0],'17':[20,21,22,23,1],'18':[20,21,22,23,2],'19':[20,21,22,23,3],'20':[16,17,18,19],'21':[16,17,18,19],'22':[16,17,18,19],'23':[16,17,18,19]}

chimera32={'0':[4,5,6,7,16],'1':[4,5,6,7,17],'2':[4,5,6,7,18],'3':[4,5,6,7,19],'4':[0,1,2,3,12],'5':[0,1,2,3,13],'6':[0,1,2,3,14],'7':[0,1,2,3,15],
'8':[12,13,14,15,24],'9':[12,13,14,15,25],'10':[12,13,14,15,26],'11':[12,13,14,15,27],'12':[8,9,10,11,4],'13':[8,9,10,11,5],'14':[8,9,10,11,6],'15':[8,9,10,11,7],
'16':[20,21,22,23,0],'17':[20,21,22,23,1],'18':[20,21,22,23,2],'19':[20,21,22,23,3],'20':[16,17,18,19,28],'21':[16,17,18,19,29],'22':[16,17,18,19,30],'23':[16,17,18,19,31],
'24':[28,29,30,31,8],'25':[28,29,30,31,9],'26':[28,29,30,31,10],'27':[28,29,30,31,11],'28':[24,25,26,27,20],'29':[24,25,26,27,21],'30':[24,25,26,27,22],'31':[24,25,26,27,23]}

#question=planted_solution(unit_cell,number_of_loops=30)
#solution=energy(question['solution'],question['ising_model'])
#solution_string=solution['string']
#solution_energy=solution['energy']

#h,J=adjgraph(question['ising_model'])

#bqm = dimod.BinaryQuadraticModel.from_ising(h, J, offset = 0.0)


#h=[]
#sampleset = dimod.ExactSolver().sample(bqm)
#sampleset.change_vartype('BINARY')
#print(sampleset.lowest())

#qpu=DWaveSampler(solver={'topology__type': 'chimera'})
#print(qpu.properties["annealing_time_range"])
#print(qpu.properties['default_annealing_time'])
#sampler = AutoEmbeddingComposite(qpu)
#sampleset = sampler.sample(bqm, num_reads=1000)
#sampleset.change_vartype('BINARY')
#print(sampleset.lowest())
#show(sampleset)

#print(solution_energy)



def short_anneal(bqm,t=False,num_reads=10):
    qpu=DWaveSampler(solver={'topology__type': 'chimera'})
    sampler = AutoEmbeddingComposite(qpu)

    #time in microseconds between 1 and 2000

    if t:
        if t< qpu.properties["annealing_time_range"][0]:
            raise ValueError('time is too short')
        elif t> qpu.properties["annealing_time_range"][1]:
            raise ValueError('time is too long')
        else:
            sampleset = sampler.sample(bqm, annealing_time=t, num_reads=num_reads)
    else:
        short_time=qpu.properties["annealing_time_range"][0]
        sampleset = sampler.sample(bqm, annealing_time=short_time, num_reads=num_reads)
    
    sampleset.change_vartype('BINARY')
    show(sampleset)

    return sampleset.lowest()

#print(short_anneal(bqm))

def reverse_annealing(bqm,reverse_schedule = [[0.0, 1.0], [5, 0.45], [99, 0.45], [100, 1.0]],initial_state='',refresh_each_run=True, num_reads=10):
    #[time, problem ratio]
    #Initial_state is a string

    if initial_state:
        initial={qubit: state for qubit,state in enumerate(list(initial_state))}
    else:
        initial={qubit: random.randint(0, 1) for qubit in list(range(len(bqm)))}
    
    
    qpu=DWaveSampler(solver={'topology__type': 'chimera'})
    sampler = AutoEmbeddingComposite(qpu)

    reverse_anneal_params = dict(anneal_schedule=reverse_schedule, initial_state=initial, reinitialize_state=refresh_each_run)
    sampleset = sampler.sample(bqm, num_reads=num_reads, **reverse_anneal_params, vartype='BINARY')
    #sampleset.change_vartype('BINARY')

    show(sampleset)
    return {'initial_state':initial_state, 'lowest_energies':sampleset.lowest()}

print(reverse_annealing(bqm))


def resacale_to_one(isingProblem):
    max_field=np.amax([abs(field) for field in isingProblem.values()])
    for coupling,strength in isingProblem.items():
        isingProblem[coupling]=strength/max_field

    return(isingProblem)



def energy_Dwave(string,couplings):

    bits=list(string)
    zeig=[1-2*int(binary) for binary in bits]
    energy=0

    for coupling,field in couplings.items():
        qubit1,qubit2=coupling
        energy=energy+zeig[qubit1]*zeig[qubit2]*field

    return {'energy':energy,'string':string}



def random_string_generator(isingProblem,trials):


    n=len(isingProblem)
    best_guess={'guess':'binary_string','energy':999}


    for number in range(trials):
        guess=format(random.randint(0,2**n-1), '0'+str(n)+'b')
        trial=energy_Dwave(guess,isingProblem)
        if trial['energy']<best_guess['energy']:
            best_guess['guess']=trial['string']
            best_guess['energy']=trial['energy']

    return {qubit:int(state) for qubit,state in enumerate(best_guess['guess'])},best_guess['energy']


def tesselationOneProblem(isingProblem, numberCopies=8, qubitNumber=24):
    dic_to_out={}

    for k in range(numberCopies):
        for coupling, field in isingProblem.items():
            qubit1,qubit2=coupling
            dic_to_out[(qubit1+k*qubitNumber,qubit2+k*qubitNumber)]=field

    return dic_to_out

def tesselationProblems(listIsingProblems,numberCopiesEach=8, qubitNumber=24):


    listUpdated=[]
    for buffer,instance in enumerate(listIsingProblems):
        holding_dic={}
        for coupling, field in instance.items():
            qubit1,qubit2= coupling
            holding_dic[(qubit1+buffer*qubitNumber,qubit2+buffer*qubitNumber)]=field
        listUpdated.append(holding_dic)

    mergedDic={instance:value for mini_dic in listUpdated for instance,value in mini_dic.items()}

    return tesselationOneProblem(mergedDic,numberCopies=numberCopiesEach,qubitNumber=qubitNumber*len(listIsingProblems))



    

raw_list_of_problems=[{(0, 4): -0.75, (0, 5): 0.5, (0, 16): 0.25,  (1, 4): -0.25,  (1, 5): -0.25,  (1, 7): 0.5, (1, 17): -0.5, (2, 4): 0.75,  (2, 5): 0.5,  (2, 6): 0.25,  (2, 7): -0.5,  (2, 18): -0.5,  (3, 5): 0.25,  (3, 19): 0.25,  (4, 12): -0.25,  (5, 13): 0.5,  (6, 14): -0.25,  (8, 12): -0.25,  (8, 13): -0.5,  (8, 14): 0.5,  (8, 15): -0.25,  (9, 12): -0.25,  (9, 14): 0.5,  (9, 15): 0.75,  (10, 12): -0.5,  (10, 13): -1.0,  (10, 14): 0.5,  (10, 15): 1.0,  (11, 12): -0.25,  (11, 14): -0.25,  (11, 15): -0.5,  (16, 20): -0.5,  (16, 21): 0.25,  (16, 23): -0.5,  (17, 20): -0.25,  (17, 21): -0.25,  (17, 22): 0.25,  (17, 23): 0.25,  (18, 20): -0.25,  (18, 23): -0.75,  (19, 21): -0.5,  (19, 22): 0.25},
{(0, 4): 1.0,  (0, 5): 2.0,  (0, 7): 2.0,  (0, 16): 1.0,  (1, 4): 1.0,  (1, 6): 1.0,  (1, 7): 4.0,  (2, 6): -2.0,  (2, 7): 1.0,  (2, 18): -1.0,  (3, 6): 3.0,  (3, 7): -1.0,  (7, 15): 2.0,  (8, 12): -4.0,  (8, 13): 3.0,  (8, 14): 4.0,  (8, 15): -3.0,  (9, 12): -1.0,  (9, 13): 2.0,  (9, 15): 1.0,  (10, 14): -3.0,  (10, 15): 3.0,  (11, 12): -1.0,  (11, 13): 1.0,  (11, 14): -1.0,  (11, 15): -1.0,  (16, 20): 1.0,  (16, 21): 2.0,  (16, 22): 1.0,  (16, 23): 1.0,  (17, 20): 1.0,  (17, 21): 1.0,  (17, 22): 2.0,  (18, 20): -2.0,  (18, 21): 1.0,  (19, 22): -1.0,  (19, 23): -1.0},
{(0, 4): -1.0,  (0, 5): -1.0,  (0, 6): -1.0,  (0, 7): 1.0,  (1, 4): -2.0,  (1, 5): 1.0,  (1, 7): -2.0,  (1, 17): 3.0,  (2, 4): -1.0,  (2, 5): 1.0,  (2, 6): 2.0,  (2, 7): -2.0,  (3, 4): -1.0,  (3, 5): -1.0,  (3, 6): 2.0,  (3, 7): -1.0,  (3, 19): -3.0,  (4, 12): 1.0,  (6, 14): 1.0,  (7, 15): 2.0,  (8, 12): 1.0,  (8, 14): -1.0,  (9, 13): -1.0,  (9, 14): -1.0,  (10, 12): -2.0,  (10, 13): 2.0,  (10, 14): 1.0,  (10, 15): -3.0,  (11, 12): -2.0,  (11, 13): 1.0,  (11, 15): -1.0,  (16, 21): -3.0,  (16, 23): -3.0,  (17, 20): -3.0,  (17, 21): 3.0,  (17, 22): 1.0,  (17, 23): 2.0,  (18, 20): 1.0,  (18, 21): -3.0,  (18, 23): -4.0,  (19, 20): 2.0,  (19, 21): -1.0,  (19, 22): -1.0,  (19, 23): 1.0},
{(0, 4): 1.0,  (0, 5): 1.0,  (0, 7): -2.0,  (1, 4): 2.0,  (1, 5): 2.0,  (1, 7): -1.0,  (1, 17): 1.0,  (2, 5): 1.0,  (2, 6): 1.0,  (2, 7): 2.0,  (3, 4): -2.0,  (3, 6): 1.0,  (3, 7): 4.0,  (3, 19): -1.0,  (4, 12): 1.0,  (7, 15): -1.0,  (8, 13): -1.0,   (8, 15): 1.0,  (9, 13): -1.0,  (9, 14): 2.0,  (9, 15): 1.0,  (10, 12): -1.0,  (10, 13): 1.0.   (11, 12): -2.0,  (11, 13): 1.0,  (11, 14): -2.0,  (11, 15): 1.0,  (16, 22): 2.0,  (16, 23): -2.0,  (17, 21): -2.0,  (17, 22): 2.0,  (17, 23): 1.0,  (18, 20): -2.0,  (18, 21): -1.0,  (18, 22): -1.0,  (18, 23): 2.0,  (19, 21): -1.0,  (19, 22): -3.0,  (19, 23): -1.0},
{(0, 4): 1.0,  (0, 5): 1.0,  (0, 6): -1.0,  (0, 7): -3.0,  (1, 4): 1.0,  (1, 5): 3.0,  (1, 6): -1.0,  (1, 7): -2.0,  (1, 17): 1.0,  (2, 5): -3.0,  (2, 6): 1.0,  (2, 7): 3.0,  (2, 18): 1.0,  (3, 4): 1.0,  (3, 5): 2.0,  (3, 6): -2.0,  (3, 7): -1.0,  (4, 12): -3.0,  (5, 13): -1.0,  (6, 14): 3.0,  (7, 15): 1.0,  (8, 12): 1.0,  (8, 15): 1.0,  (9, 12): -2.0,  (9, 13): -2.0,  (9, 14): 1.0,  (9, 15): 3.0,  (10, 12): 2.0,  (10, 13): -1.0,  (10, 15): 1.0,  (11, 12): -2.0,  (11, 15): -2.0,  (16, 20): -1.0,  (16, 21): 2.0,  (16, 23): -1.0,  (17, 20): -4.0,  (17, 22): -3.0,  (17, 23): 2.0,  (18, 20): 1.0,  (18, 21): -1.0,  (18, 23): -1.0,  (19, 20): -2.0,  (19, 21): 3.0,  (19, 22): 1.0},
{(0, 4): 5.0,  (0, 5): 1.0,  (0, 6): 1.0,  (0, 7): -2.0,  (0, 16): -1.0,  (1, 5): -1.0,  (1, 6): -1.0,  (1, 7): 3.0,  (1, 17): -1.0,  (2, 5): -1.0,  (2, 6): -1.0,  (2, 7): -1.0,  (2, 18): 1.0,  (3, 4): 3.0,  (3, 6): -2.0,  (3, 19): 1.0,  (4, 12): -2.0,  (5, 13): 1.0,  (6, 14): 1.0,  (7, 15): -4.0,  (8, 12): 1.0,  (8, 13): 1.0,  (8, 15): -2.0,  (9, 12): -1.0,  (9, 13): 3.0,  (9, 14): -1.0,  (9, 15): 3.0,  (10, 12): 4.0,  (10, 13): -1.0,  (10, 14): 1.0,  (10, 15): -2.0,  (11, 12): -2.0,  (11, 14): -1.0,  (11, 15): -1.0,  (16, 20): 1.0,  (16, 22): -2.0,  (17, 21): -2.0,  (17, 22): 1.0,  (18, 20): -2.0,  (18, 22): 2.0,  (18, 23): -1.0,  (19, 20): 1.0,  (19, 22): 1.0,  (19, 23): -1.0},
{(0, 4): 2.0,  (0, 5): -1.0,  (0, 16): -1.0,  (1, 4): 1.0,  (1, 5): -2.0,  (1, 6): 3.0,  (1, 7): 2.0,  (2, 6): 2.0,  (2, 7): 1.0,  (2, 18): 1.0,  (3, 4): -1.0,  (3, 5): -2.0,  (3, 7): 1.0,  (5, 13): 1.0,  (6, 14): 1.0,  (8, 12): 1.0,  (8, 13): 1.0,  (9, 12): -2.0,  (9, 14): 2.0,  (9, 15): -2.0,  (10, 13): -1.0,  (10, 14): 1.0,  (10, 15): -4.0,  (11, 12): -1.0,  (11, 13): -1.0,  (11, 14): 2.0,  (16, 21): 1.0,  (16, 22): 1.0,  (16, 23): -3.0,  (17, 20): 1.0,  (17, 22): -1.0,  (17, 23): 2.0,  (18, 21): -1.0,  (19, 20): -1.0,  (19, 21): 4.0,  (19, 22): 4.0,  (19, 23): -1.0},
{(0, 6): -1.0,  (0, 16): 1.0,  (1, 4): -1.0,  (1, 5): -1.0,  (1, 6): 1.0,  (1, 17): -1.0,  (2, 5): -1.0,  (2, 6): 1.0,  (2, 7): 2.0,  (3, 4): 3.0,  (3, 5): -1.0,  (3, 6): 2.0,  (3, 7): 4.0,  (5, 13): -1.0,  (6, 14): 1.0,  (7, 15): -2.0,  (8, 12): -3.0,  (8, 13): 1.0,  (8, 15): -2.0,  (9, 12): -1.0,  (9, 13): 2.0,  (9, 14): 1.0,  (10, 12): 1.0,  (10, 13): 1.0,  (10, 15): -2.0,  (11, 12): -1.0,  (11, 13): 1.0,  (11, 15): -2.0,  (16, 20): 2.0,  (16, 21): 2.0,  (16, 22): 1.0,  (17, 20): 2.0,  (17, 21): -1.0,  (17, 22): -1.0,  (17, 23): -3.0,  (18, 20): 2.0,  (18, 21): -1.0,  (18, 23): -1.0,  (19, 21): 2.0,  (19, 22): -2.0,  (19, 23): -2.0},
{(0, 5): -1.0,  (0, 6): -3.0,  (0, 7): -1.0,  (0, 16): 1.0,  (1, 4): 3.0,  (1, 5): 2.0,  (1, 6): 3.0,  (1, 7): 1.0,  (1, 17): 1.0,  (2, 4): 2.0,  (2, 6): 2.0,  (2, 7): -1.0,  (2, 18): 1.0,  (3, 4): -1.0,  (3, 5): 2.0,  (3, 6): 2.0,  (3, 19): -1.0,  (5, 13): -1.0,  (7, 15): 1.0,  (8, 12): 2.0,  (8, 13): -5.0,  (8, 14): 4.0,  (8, 15): 1.0,  (9, 12): -2.0,  (9, 13): 3.0,  (9, 15): -1.0,  (10, 12): 1.0,  (10, 15): 1.0,  (11, 12): -1.0,  (11, 13): 1.0,  (16, 20): -2.0,  (16, 21): -3.0,  (16, 22): 2.0,  (16, 23): -2.0,  (17, 20): 1.0,  (17, 21): 1.0,  (17, 23): -1.0,  (18, 20): 1.0,  (18, 22): 1.0,  (18, 23): -1.0,  (19, 20): -4.0,  (19, 22): 1.0,  (19, 23): 2.0},
{(0, 4): -2.0,  (0, 5): 2.0,  (0, 6): 1.0,  (0, 7): 1.0,  (1, 4): -1.0,  (1, 5): -1.0,  (1, 6): -2.0,  (1, 7): -2.0,  (2, 4): 1.0,  (2, 5): -1.0,  (2, 7): -1.0,  (2, 18): 1.0,  (3, 4): -3.0,  (3, 5): 4.0,  (3, 6): 1.0,  (3, 7): 1.0,  (3, 19): 1.0,  (4, 12): -1.0,  (7, 15): 1.0,  (8, 12): 1.0,  (8, 13): -1.0,  (8, 14): -2.0,  (9, 12): 3.0,  (9, 13): -4.0,  (9, 14): 1.0,  (9, 15): 2.0,  (10, 12): 1.0,  (10, 13): 4.0,  (10, 15): -1.0,  (11, 13): -1.0,  (11, 14): -1.0,  (16, 20): 4.0,  (16, 22): -2.0,  (17, 20): 1.0,  (17, 21): -1.0,  (18, 21): 2.0,  (18, 22): 1.0,  (18, 23): -2.0,  (19, 20): 3.0,  (19, 21): -1.0,  (19, 22): -3.0}]

list_of_problems=[resacale_to_one(problem) for problem in raw_list_of_problems]

problem_to_run=tesselationProblems(list_of_problems)






    
        











