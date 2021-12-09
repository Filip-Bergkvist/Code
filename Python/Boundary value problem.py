import matplotlib.pyplot as plt
import numpy as np

# Nonlinear BVP from lab assignment
def f(x,y):
    '''
    Second order differential equation equals return-object
    '''
    return (3 / 2) * y

# Functions jacobi, runge_kutta and Secant are used from my work on earlier labs

def runge_kutta(x0,y0, h, end_interval, function,alpha):
    '''
    Parameters
    ----------
    x0 : start interval x-coordinate
    y0 : starting y-value
    h : step-size 
    end_interval : ending x-coordinate interval
    function: differential equation
    alpha: function(x0) (lower boundary value)

    Returns
    -------
    result : list of interpolated(x,y) values using Runge-Kutta
    '''
    
    result = []
    
    n = (end_interval - x0) / h 
    
    result.append((x0,alpha))
    x0 = round(x0 + h, 5)

    for i in range(int(n)):
        
        k_1 = h * function(x0,y0)
        k_2 = h * function(x0+(h/2) , y0 + (k_1/2))
        k_3 = h * function(x0 + (h/2), y0 + (k_2/2))
        k_4 = h * function((x0 + h), (y0 +k_3))
                    
        y = y0 + (k_1 / 6 + k_2 / 3 + k_3 / 3 + k_4 / 6)
        
        result.append((x0,y))
        
        x0 = round(x0 + h, 5)
        y0 = y

    #print("The approximation of " + "y("  + str(end_interval) + ") with Runge-Kutta, " + "with h = " + str(h) + " is: " + str(result[int(n)][1]) + "\n")
    
    return result

def Secant(a, b, e, imax, f ):
    ''' 
    a = first point interval
    b = second point interval
    e = tolerance
    imax = max iterations
    f = function 
    '''
    k = 1
    
    value_b = f(b)
    value_a = f(a)
    
    while k < imax: 
        
        root = (a * f(b) - b * f(a)) / (f(b)- f(a))
        
        value = f(root)
        
        
        if abs(value) < e:
            break
        
        elif abs(b-a) < e:
            break
        
        elif value_b * value > 0:
            b = root
            value_b = value
            
        elif value_a * value > 0:
            a = root
            value_a = value
 
        k = k + 1
        
    return root

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
        
    #print("Solution to linear equation using jacobi: \n")
    #for i in range(len(b)):
     #   print("X%d = %0.3f" % (i+1,x[i]))
    
    return x

#Shooting method 

def shooting(x0,function, delta_x, end_interval, max_iterations, alpha, beta):
    '''
    
    Parameters
    ----------
    x0 : lower boundary 
    function : second order differential equation
    delta_x : step size
    end_interval : upper boundary
    max_iterations : max number of iterations
    alpha : lower boundary value
    beta : upper boundary value

    Returns
    -------
    Approximation of x(t) value with shooting method
    using Runge-Kutta and Secant method.

    '''
    
    #Initial guess
    x_prim = (beta - alpha) / (end_interval - x0)
    
    table = []
    
    iterations = 0
    
    # Function to return the value x(upper boundary)
    def uppervalue_runge(x0, x_prim, delta_x, end_interval ,function,alpha):
        
        lista = runge_kutta(x0, x_prim, delta_x, end_interval ,function,alpha)
        
        return lista[-1][1]
    
    # This function is used by the Secant method in the 
    # iteration below to find root and get a better "x_prim" guess. 
    def helper_rungekutta(x_prim):
        
        return runge_kutta(x0, x_prim, delta_x, end_interval ,function,alpha)[-1][1] - beta
    
    
    for i in range(max_iterations):
        
        result = uppervalue_runge(x0, x_prim, delta_x, end_interval ,function,alpha)

        iterations = iterations+ 1
        
        # Save all the Runge-Kutta approximations in a table to plot later
        curve = runge_kutta(x0, x_prim, delta_x, end_interval ,function,alpha) 
        
        table.append((curve))
        
        if abs(result - beta) < 0.0001:
            
            break
        
        else: 
            
           # Using Secant to solve root to function result - beta = 0, create new guess for 
           # x'(0). "result" is a function of x_prim (value received with help of runge_kutta)
           x_prim = Secant(0, 1, 0.01, 10, helper_rungekutta)
           
    # Iteration below to plot all lines of runge-kutta method above. (All "shots")
    k = 1 
    for i in table:
        plt.title("Shooting method approximation with delta x = " + str(delta_x))
        plt.ylabel("x(t)")
        plt.xlabel("t" )
        plt.plot([x[0] for x in i], [x[1] for x in i] , label = "Attempt " + str(k))
        plt.legend(loc= "upper right")
        
        k = k + 1
    
    print("The approximation of " + "x("  + str(end_interval) + ") with Shooting method using Runge-Kutta, " + "with delta x = " + str(delta_x) + ", number of iterations =  "+ str(iterations) + ", and final guess x'(0) = " + str(x_prim) + " is: " + str(result)) 
    
    return "\n" 

print(shooting(0, f, 0.1, 1, 10,4,1))
plt.show()
print(shooting(0, f, 0.01, 1, 10,4,1))
plt.show()

# Finite difference method 

def finite_difference(t0, end_interval, delta_x, alpha, beta, function):
    '''
    

    Parameters
    ----------
    t0 : lower boundary
    end_interval : upper boundary
    delta_x : step size
    alpha : lower boundary value
    beta : upper boundary value
    function : second order differential equation in the lab assignment

    Returns
    -------
    We split up the BVP in to matrix form A X = B and use defined jacobi method 
    to solve this sytem of linear equations. Graph of result.

    '''

    n = int((end_interval - t0) / delta_x)
    
    # Initialize matrix A 
    A = np.zeros((n+1,n+1))
    A[0,0] = 1
    A[n,n] = 1
    
  
    #Use the definition of second order differential equation
    for i in range(1,n):
        A[i,i-1] = 1
        A[i,i] = -2 
        A[i,i+1] = 1
    
    # Here we adjust the matrix A according to the equation in the lab assignment
    # This would not be correct if the differential equation problem was different
    x_t = function(1,1)
    for i in range(1,n):
        A[i,i] = A[i,i] -  x_t * delta_x ** 2 * (3/2)
    
    # Initialize matrix B. Boundary values and rest zeros. 
    #(Also depending on the differential eq.)
    B = np.zeros(n+1)
    B[0] = alpha
    B[n] = beta
    
    for i in range(1,n):
        B[i] = 0
    
    # Solving the linear equation using jacobi
    result = jacobi(A,B)
    
    # All t values
    t = []
    for i in range(n + 1):
        t.append(round(i * delta_x,3))
        
    plt.title("Finite difference method")
    plt.xlabel("t" )
    plt.ylabel("x(t)")
    plt.plot(t,result, label = "Delta x: " + str(delta_x))
    plt.legend(loc = "upper right")
    
    return ""


print("Check graph for result of finite difference method! ")
print(finite_difference(0, 1, 0.1, 4, 1,f))
print(finite_difference(0, 1, 0.01, 4, 1,f))
plt.show()
