import numpy as np

A = np.array([
    [0,1,1],
    [1,1,2],
    [0,0,3]])

def QR(a): 
    
    '''
    Parameters
    ----------
    a : n x n matrix

    Returns
    -------
    QR decomposition for matrix R based on 
    Gram-Schmidt method. 

    '''
    
    # All zeros in R matrix at first
    R = np.copy(a)
    for i in range(len(a)):
        for j in range(len(a)):
            R[i][j] = 0
    
    # Q identical to input a at first 
    Q = np.copy(a)

    
    for i in range(len(a)):
        
        # Norm vector
        summa = 0
        for k in range(len(a)):
            
            summa += a[k,i] ** 2
            
        R[i,i] = summa ** (1/2)
            
    
        Q[:,i] = Q[:,i] / R[i,i]
        

        for j in range(i+1, len(a)):
            
            # Dot product 
            summa_dot = 0
            for k in range(len(a)):
                
                summa_dot += Q[k,i] * a[k,j]
                
            R[i,j] = summa_dot 
            
            Q[:,j] = Q[:,j] - R[i,j] * Q[:,i]
            
    print("Input Matrix A: \n" + str(a))
    print("Matrix Q: \n" + str(Q))
    print("Matrix R: \n" + str(R))
    
    
    return "Done"
    

print(QR(A))
