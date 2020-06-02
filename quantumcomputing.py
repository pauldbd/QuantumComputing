import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

        
class QuantumCom_Vector:

    def __init__(self):
        #The matrices for the functions
        #Pauli matrices
        self.sigmax = [[0,1], [1,0]]
        self.sigmay = [[0, -1j], [1j, 0]] # J --> i = sqrt(-1)
        self.sigmaz = [[1,0],[0,-1]]

        #Controlled Gates
        self.controlledx1 = [[1,0,0,0], [0,1,0,0], [0,0,0,1], [0,0,1,0]]
        self.controlledx2 = [[1,0,0,0], [0,0,0,1], [0,0,1,0], [0,1,0,0]]
        #General Functions
        sqr2 = np.sqrt(2)
        self.hadamard = [[1/sqr2, 1/sqr2], [1/sqr2, -1/sqr2]]

        self.functions = []
        self.inputs = []

    #The following functions will flip the arbitrary vector 180 degrees in certain axis which is the name of the functions
    def x(self, qbit):
        qbit.vector = list(np.dot(qbit.vector, self.sigmax))

    def y(self, qbit):
        qbit.vector = list(np.dot(qbit.vector, self.sigmay))

    def z(self, qbit):
        qbit.vector = list(np.dot(qbit.vector, self.sigmaz))

    #General Functions such as the Hadamard Gate
    def Hadamard(self, qbit):
        qbit.vector = list(np.dot(qbit.vector, self.hadamard))
    #Controlled gates
    def cx(self, q):
        x = q[0].vector
        y = q[1].vector

        v = self.tensorVector(x, y)
        return list(np.dot(v, self.controlledx))

    #Measuring the qbit -> Getting the probability density of the wavefunction
    def measure(self, qbit):
        #P = []
        #P.append([np.abs(qbit.vector[0])**2, np.abs(qbit.vector[1]**2)])

        P = []
        for i in range(len(qbit.vector)):
            P.append(round(np.abs(qbit.vector[i])**2, 2))

        return P
        
    def addCircuit(self, functions, inputs):
        self.functions.append(functions)
        self.inputs.append(inputs)

    def runCircuit(self, returnLast=False):
        # Every function is made such that there would only have to be one input: if there is a need for a multiple inputs
        # you will be able to provide them in an array
        if returnLast == True:
            for i in range(len(self.functions) - 1):
                self.functions[i](self.inputs[i])

            return self.functions[len(self.functions) - 1](self.inputs[len(self.functions) - 1])
        else:
            for i in range(len(self.functions)):
                self.functions[i](self.inputs[i])

    def TensorProductVector(self, x, y):
        #Tensor product of two, two dimensional vectors
        #Create a vector with 2^(size of x) -> 4 
        
        if type(x) != int and type(y) != int:
            vector = list(np.zeros(np.size(x) * np.size(y)))
            d = 0
            for i in range(len(x)):
                for j in range(len(y)):
                    vector[d] = (x[i]*y[j])
                    d += 1
            
            return vector

        else:
            return x*y

    def BellState(self, qubits):
        a = qubits[0]
        b = qubits[1]

        self.Hadamard(a)
        return self.cx([a, b])

    def entangledGate(self, qubit, gate, qubitnum, returnMatrix=False):
        m = 1
        identityFactor = len(qubit)/np.log2(np.size(gate))/(np.log2(len(qubit))-np.log2(np.shape(gate)[0]))
        for i in range(len(qubit)):
            if i == qubitnum:
                m = np.kron(m, gate)

            else:
                m = np.kron(m, np.identity(identityFactor))
        if returnMatrix:
            return np.dot(m, qubit),m

        return np.dot(m, qubit)

        

    def normalize(self, qubit):
        a = qubit.vector[0]
        b = qubit.vector[1]

        S = np.abs(a) ** 2 + np.abs(b) ** 2
        m = np.sqrt(1/S)

        qubit.vector = list(np.dot(qubit.vector, m))

    def showBlochRepresentation(self, qubit):
        #qubit is vector here
        a = qubit[0]
        b = qubit[1]

        theta = 2*np.arccos(a)
        x = b/(np.sin(theta/2))

        phi = np.arcsin(x.imag)
        return theta, phi
        
class qubit:
    def __init__(self, vector=[1,0]):
        self.vector = vector

    def group(self, *qubits):
        tensorPrd = 1
        for qubit in qubits:
            tensorPrd = QuantumCom_Vector().TensorProductVector(tensorPrd, qubit.vector) 

        return tensorPrd

    def createQubit(self, q):

        res = q[::-1] 

        for i in range(len(q)):
            if q[i] == 0:
                q[i] = [1,0]

            else:
                q[i] = [0,1]


        d = 1
        for i in range(len(q)):
            d = np.kron(d, q[i])

        return d

    
        
class visualizeQubit:
    def __init__(self):
        self.BlochX = None


    def Bloch2D(self, qubitVector):
        #Hello
        fig = plt.figure()
        ax = plt.axes()

        plt.axis("equal")


        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')

        x = np.cos(np.arange(0, 10, 0.00001))
        y = np.sin(np.arange(0, 10, 0.00001))
        ax.plot(x, y, color='k')

        ax.scatter(0,0, s=25, color='k')

        ax.arrow(0,0,0,1, width=0.01,length_includes_head=True, head_length=0.1, color='k')
        ax.arrow(0,0,1,0, width=0.01,length_includes_head=True, head_length=0.1, color='k')
        ax.arrow(0,0,0,-1, width=0.01,length_includes_head=True, head_length=0.1, color='k')
        ax.arrow(0,0,-1,0, width=0.01,length_includes_head=True, head_length=0.1, color='k')
        plt.text(0.05,.85, '|1>', fontsize=12)
        plt.text(.76,0.07, '|0>', fontsize=12)

        x = []
        y = []
        for i in range(len(qubitVector)):
            x.append(qubitVector[0])
            y.append(qubitVector[1])

        ax.arrow(0,0,x,y, width=0.015,length_includes_head=True, head_length=0.1, color='k')

        X = []
        Y = []

        for i in range(len(x)):
            if x[i] >= 0:
                X.append(x[i]+0.08)
            else:
                X.append(x[i]-0.23)

        for i in range(len(y)):
            if y[i]>=0:
                Y.append(y[i]+0.03)

            else:
                Y.append(y[i]-0.05)

        plt.text(X,Y, '|Ψ>', fontsize=12)

        plt.show()

    def Bloch3D(self, qubitVectors, showTransformation=False, colors=[]):
        
        if colors == []:
            for i in range(len(qubitVectors)):
                colors.append('k')
            

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("equal")

        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        ax.plot_surface(x, y, z, color="k", alpha=0.2)
        ax.plot_wireframe(x, y, z, color="k", alpha=0.5)

        
        a = np.arange(-1, 1, 0.01)

        n = len(a)
        b = np.zeros(n)

        ax.plot(a, b, b, color='k', linewidth=2.5, alpha=0.7)
        ax.plot(b, a, b, color='k', linewidth=2.5, alpha=0.7)
        ax.plot(b, b, a, color='k', linewidth=2.5, alpha=0.7)

        ax.set_xlabel('|+> and |-> Basis(X)')
        ax.set_ylabel('Imaginary Numbers(Y)')
        ax.set_zlabel('Standard Basis(Z)')

        ax.text(0,0,1.1, "|0>", fontsize=13)
        ax.text(0,0,-1.2, "|1>", fontsize=13)

        ax.text(1.1, 0,0, "|+>", fontsize=13)
        ax.text(-1.3, 0,0, "|->", fontsize=13)

        ax.text(0, 1.1,0, "i", fontsize=13)
        ax.text(0, -1.2,0, "-i", fontsize=13)

        '''
        Note to myself: 
        The angles such as theta can be represented by taking the cosine of the the position on z axis
        This could possibily be used for the phi as well!
        '''

        X = []
        Y = []
        Z = []

    
        d = 0
        for qubitVector in qubitVectors:
            theta, phi = QuantumCom_Vector().showBlochRepresentation(qubitVector)
            
            x = np.cos(phi) * np.sin(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(theta)

            X.append(x)
            Y.append(y)
            Z.append(z)

            ax.quiver(0,0,0,x, y, z, color=colors[d], linewidths=3)
            d +=1 

        ax.plot(X, Y, Z, '--', linewidth=3, color='k')
        

        plt.show()

class QuantumCom_BraKet:
    def __init__(self):
        pass
    def VectorToKet(self, vector):
        n = np.size(vector)
        b = {} 
        
        for i in range(n):
            if np.log2(n) - len(bin(i)[2:]) == 0:
                b[bin(i)[2:]]=vector[i]
            else:
                
                a = np.log2(n) - len(bin(i)[2:])
                d = ''
                for c in range(int(a)):
                    d+='0'
                b[d + bin(i)[2:]]=vector[i]
                
                

        return b
 
    def KetToVector(self, braket):
        d = list(np.zeros(len(braket)))
        for keys in braket:
            b = 1
            for i in keys:
                if i == 1:
                    a=[0,1]
                else:
                    a = [1,0]
                b = np.kron(b, a)
            
            #d += b*braket[keys]
            d = np.add(d, b*braket[keys])
        return list(d)

qc = QuantumCom_Vector()
q1 = qubit([0,1])
q2 = qubit([1,0])
q3 = qubit([1,0])

q1_2 = qubit()
q1_2.vector = qc.BellState([q1, q2])

QubitSystem = qubit(np.kron(q1_2.vector, q3))

qc.entangledGate(QubitSystem, )


