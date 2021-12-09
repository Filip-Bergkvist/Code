mport numpy as np



B = np.array([[-5.173],
      [-5.458],
      [4.415]])

A = np.array([[0.143,0.357,2.01],
      [-1.31,0.911,1.99],
      [11.2,-4.3,-0.605]])

        

def Gaussian_pivoting(a, b):
    '''
    

    Parameters
    ----------
    a : n x n matrix
    b : n x 1 vector

    Returns
    -------
    Solution to the linear equation ax = b. 
    Uses Gaussian elimination with pivoting. 
    Each step in the elimination with pivoting 
    is printed out.

    '''
    
    print("Gaussian elimination with pivoting: \n")
    
    print("Input matrix left: " + str(a) + "\n")
    print("Input matrix right: " + str(b) + "\n")
    counter = 1
    for k in range(len(a)-1):
        
        for i in range(k+1, len(a)):
            
            # Pivoting takes place
            if abs(a[i,k]) > abs(a[k,k]):
                
                a[[k,i]] = a[[i,k]]
                
                b[[k,i]] = b[[i,k]]
                          
                
        for i in range(k+1, len(a)):
            if a[i][k] == 0: 
                continue
                
            m = a[k][k] / a[i][k]
                
            
            for j in range(k, len(a)):
                    
                a[i][j] = a[k][j] - a[i][j]*m
                
            b[i][0] = b[k][0] - b[i][0] * m
            
            print("Step " + str(counter) + " " + "left side : " + str(a) + "\n")
            print("Step " + str(counter) + " " + "right side : " + str(b) + "\n")
            counter += 1
      
    b[(len(a)-1)][0] = b[(len(a)-1)][0] / a[len(a)-1][(len(a)-1)]
    
    
    for i in range(len(a)-2 ,-1, -1):
        
        for j in range(i+1, len(a)):
            
            b[i][0] = b[i][0] - a[i][j]*b[j][0] 
            
        b[i][0] = b[i][0] / a[i][i]
    
    print("Answer to equation: \n")
    for i in range(len(b)):
        print("X%d = %0.5f" % (i+1,b[i][0]))
    
    return "Done"



print(Gaussian_pivoting(A, B))
