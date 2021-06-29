import numpy as np
import random

from numpy.lib.twodim_base import tri

def adjgraph(adj):
    adj=np.array(adj)
    uptri=np.triu(adj,1)
    
    h={int(qubit):singlefield for qubit,singlefield in enumerate(np.diagonal(adj))}
    j={(qubit1,qubit2):interaction for qubit1,conn in enumerate(uptri) for qubit2,interaction in enumerate(conn) if interaction!=0}

    return h,j

mat=[[1,2,3],[4,5,6],[7,8,9]]
#print(adjgraph(mat))

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


#unit_cell={'0':[4,5,6,7],'1':[4,5,6,7],'2':[4,5,6,7],'3':[4,5,6,7],'4':[0,1,2,3],'5':[0,1,2,3],'6':[0,1,2,3],'7':[0,1,2,3]}

#question=planted_solution(unit_cell,number_of_loops=6,scale=False)
#print(question['ising_model'])
#print(adjgraph(question['ising_model']))


def energy(string, adj):
    n=len(adj)
    h=np.diag(adj)
    adj=np.array(adj/2)-np.diag(h)/2

    bits=list(string)
    zeig=np.array([1-2*int(x) for x in bits])
    value=np.dot(zeig,np.dot(adj,zeig))+np.dot(np.array(h),zeig)
    

    return {'energy':value,'string':string}



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

    print(mergedDic)

    return tesselationOneProblem(mergedDic,numberCopies=numberCopiesEach,qubitNumber=qubitNumber*len(listIsingProblems))

def relabelChimera(connection_dictionary):
    
    remapped_dictionary={}
    for connection,field in connection_dictionary.items():
        qubit1,qubit2=connection
        if 0<= qubit1 <= 3:
            qubit1=qubit1+4
        elif 4<= qubit1 <= 7:
            qubit1=qubit1-4
        elif 16<= qubit1 <= 19:
            qubit1=qubit1+4
        elif 20<= qubit1 <= 23:
            qubit1=qubit1-4

        if 0<= qubit2 <= 3:
            qubit2=qubit2+4
        elif 4<= qubit2 <= 7:
            qubit2=qubit2-4
        elif 16<= qubit2 <= 19:
            qubit2=qubit2+4
        elif 20<= qubit2 <= 23:
            qubit2=qubit2-4

        remapped_dictionary[(qubit1,qubit2)]=field

    return remapped_dictionary
        
