import numpy as np

def adjgraph(adj):
    adj=np.array(adj)
    uptri=np.triu(adj,1)
    
    h={int(qubit):singlefield for qubit,singlefield in enumerate(np.diagonal(adj))}
    j={(qubit1,qubit2):interaction for qubit1,conn in enumerate(uptri) for qubit2,interaction in enumerate(conn) if interaction!=0}

    return h,j

mat=[[1,2,3],[4,5,6],[7,8,9]]
print(adjgraph(mat))

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


unit_cell={'0':[4,5,6,7],'1':[4,5,6,7],'2':[4,5,6,7],'3':[4,5,6,7],'4':[0,1,2,3],'5':[0,1,2,3],'6':[0,1,2,3],'7':[0,1,2,3]}

question=planted_solution(unit_cell,number_of_loops=6,scale=False)
print(question['ising_model'])
print(adjgraph(question['ising_model']))
