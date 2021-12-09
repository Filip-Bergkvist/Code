A_1 = [[1,5],
     [-2,-7]]


B_1 = [[7],
     [-5]]

A_2 = [[1,2,4],
      [3,8,14],
      [2,6,13]]

B_2 = [[3],
      [13],
      [4]]


def Gaussian(a, b):  
    '''
    Parameters
    ----------
    a : n x n matrix
    b : n x 1 vector
        
    Returns
    -------
    Solution to the linear equation ax = b. 
    Uses Gaussian elimination.
        
    '''
    
    
    print("Gaussian for linear eq :\n" + "Left:" +  str(a) + "\n" + "Right:" +  str(b))
    m = 0

    for i in range(len(a)):
        
        for j in range(i+1, len(a)):
            
            try:
                m = a[j][i] / a[i][i]
            
            except ZeroDivisionError:
                print("Division by zero, answer not to be trusted")
                break
            
            b[j][0] = b[j][0] - (m*b[i][0]) 
            for k in range(len(a)):
                a[j][k] = a[j][k] - m*a[i][k]
                
    
    # Get the value of bottom right variable in matrix a
    # and assigning it to b, b will represent the variable-values in the 
    # solution  
    b[(len(a)-1)][0] = b[(len(a)-1)][0] / a[len(a)-1][(len(a)-1)]
    
    try: 
        for i in range(len(a)-2 ,-1, -1):
    
            
            for j in range(i+1, len(a)):
    
                
                b[i][0] = b[i][0] - a[i][j]*b[j][0] 
                
            
            
            b[i][0] = b[i][0] / a[i][i]
            
    except ZeroDivisionError:
        print("Division by zero, answer not to be trusted")
        
        
    for i in range(len(b)):
        print("X%d = %0.2f" % (i+1,b[i][0]))
    
    return "Done \n"



print(Gaussian(A_1, B_1))

print(Gaussian(A_2,B_2))
