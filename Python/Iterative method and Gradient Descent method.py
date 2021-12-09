import numpy as np

A = np.array([[5,-2,3],
             [-3,9,1],
             [2,-1,-7]])


B = np.array([-1,2,3])

def jacobi(a, b):
    '''

    Parameters
    ----------
    a : n x n matrix
    b : vector

    Returns
    -------
    A solution x to the linear equations 
    ax = b

    '''

    D= []
    for i in range(len(a)):
        D.append(a[i,i])
        a[i,i] = 0
        
    x = np.zeros(len(a[0]))

    
    for k in range(30):
       x = (b - np.dot(a,x)) / D
        
    
    print("Solution to linear equation using jacobi: \n")
    for i in range(len(b)):
        print("X%d = %0.3f" % (i+1,x[i]))
    
    return "Done \n"

print(jacobi(A,B))


''' 

Using function and gradient defined in the lab assignment below. With values
on a = 1 and b = 5. Also using the start [-1.4,2] as used in
the lab assignment.

'''
def f(x):
    
    return (1-x[0]) ** 2 + 5 * (x[1] - (x[0]**2))**2

def gradient(x):
    
    return np.array([-20 * x[0] * (x[1]-(x[0]**2)) - 2 * (1 - x[0]), 
                     10 * (x[1] - (x[0]**2))])


def gradient_descent(function,gradient, e, M, gamma):
    '''

    Parameters
    ----------
    function : function to minimize
    
    gradient : gradient of above function
        
    e : Tolerance(precision)
    
    M : Max number of iterations
    
    gamma : learning step size

    Returns
    -------
    The minimum value and the minimum point of function with
    a certain starting value.

    '''
    
    # Define start points
    start_value = np.array([-1.4, 2])
    
    i = 0
    
    while i < M:
        
        vector = -gamma * gradient(start_value)
        
        if abs(vector[0]) < e: 
            if abs(vector[1]) < e:
                break
        
        start_value += vector
        
        i = i + 1
    
    
    
    print("The minimum value is: \n" + str(function(start_value)) + "\n")
    print("The minimum point is: \n" + str(start_value) + "\n")
    
    return "Done"


print(gradient_descent(f, gradient, 0.00001, 1000, 0.03))
