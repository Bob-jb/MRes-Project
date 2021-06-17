import numpy as np
import dimod

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

question=planted_solution(unit_cell,number_of_loops=20,scale=False)
solution=energy(question['solution'],question['ising_model'])
solution_string=solution['string']
solution_energy=solution['energy']

h,J=adjgraph(question['ising_model'])

bqm = dimod.BinaryQuadraticModel.from_ising(h, J, offset = 0.0)
h=[]
sampleset = dimod.ExactSolver().sample(bqm)
sampleset.change_vartype('BINARY')
print(sampleset.lowest())












