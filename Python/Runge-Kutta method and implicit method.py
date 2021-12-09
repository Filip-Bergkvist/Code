import matplotlib.pyplot as plt
import math

# Assignment 1

def f(t,y):
    
    return (1 / (t ** 2)) - (y / t) - y ** 2

def exact_solution(x): 

    return  -1 / x

def euler_method(x0,y0,h):
    '''
    
    Parameters
    ----------
    x0 : starting x-value
    y0 : starting y-value
    h : step-size

    Returns
    -------
    result : list of tuples (x,y) with interpolated values using Euler method. 

    '''
    
    result = []
    
    n = 1 / h
    
    y_prim = 0
    
    for i in range(int(n+1)):
        
        y_prim = f(x0,y0)
        
        y = y0 + h * y_prim
        
        result.append((x0,y))
        
        x0 = round(x0 + h, 5)
    
        y0 = y
    
    print("The approximation of y(2) " + "with h = " + str(h) + " is: " + str(result[int(n)][1]) + "\n")
    
    return result

print("Assignment 1a): \n")

data_a = euler_method(1,-1,0.05)

print("Exact solution: " + str(exact_solution(2)) + "\n")

Query_points = [1.052, 1.555, 1.978]

def linear_splines(data,x_q,metod):
    '''

    Parameters
    ----------
    data : list of tuples with (x,y) values 
    x_q : query points (x coordinates)
    metod : string of which method is being used, for label in plot. 
    Returns
    -------
    Interpolated values for the query points
    
    '''
    
    result = []
    
    for i in x_q:

        for k in range(len(data)):
            
            output = 0
            
            if i == data[k][0]:
                
                result.append((i, data[k][1]))
            
            if i > data[k][0] and i < data[k+1][0]:
               
                output = data[k][1] + (i - data[k][0]) * ((data[k+1][1] - data[k][1]) 
                                                        / (data[k+1][0] - data[k][0]))
                result.append((i, output))
                
                break
    
    print("Interpolated values linear: \n")
    for i in range(len(result)):
        print("X coordinate: " + " " + str(result[i][0]) + " " + "Interpolated value: " + str(result[i][1])+ "\n")
    
    
    x_val = [x[0] for x in result]
    y_val = [x[1] for x in result]
    
    plt.title("Linear " + str(metod))
    plt.plot(x_val, y_val)
    plt.show()
    return "Done \n"

print("Assignment 1b): \n")
print(linear_splines(data_a, Query_points, "Euler"))

print("Exact solutions: \n")
print("1.052: " + str(exact_solution(1.052)))
print("1.555: " + str(exact_solution(1.555)))
print("1.978: " + str(exact_solution(1.978)) + "\n")


def euler_modified(x0,y0,h):
    '''
    Parameters
    ----------
    x0 : starting x-value
    y0 : starting y-value
    h : step-size 

    Returns
    -------
    result : list of tuples (x,y) with interpolated values using modified Euler method

    '''
    
    result = []
    
    n = 1 / h
    
    y_prim = 0
    
    for i in range(int(n+1)):
        
        y = y0 + h * f(x0,y0)
        
        y_prim = y0 + h/2 * (f(x0,y0) + f(x0 + h,y))
        
        result.append((x0,y_prim))
        
        x0 = round(x0 + h, 5)
    
        y0 = y
    
    print("The approximation of y(2) euler modified " + "with h = " + str(h) + " is: " + str(result[int(n)][1]) + "\n")
    
    return result

print("Assigment 1c): \n")

data_c = euler_modified(1, -1, 0.05)
print("Exact solution: " + str(exact_solution(2)) + "\n")

print("Assignment 1d): \n")

print(linear_splines(data_c, Query_points, "Euler modified"))

print("Exact solutions: \n")
print("1.052: " + str(exact_solution(1.052)))
print("1.555: " + str(exact_solution(1.555)))
print("1.978: " + str(exact_solution(1.978)) + "\n")


print("Assignment 1e): \n")

def runge_kutta(x0,y0, h, end_interval, function):
    '''
    
    Parameters
    ----------
    x0 : start interval x-coordinate
    y0 : starting y-value
    h : step-size 
    end_interval : ending x-coordinate interval

    Returns
    -------
    result : list of tuples(x,y) with interpolated values using Runge-Kutta

    '''
    
    result = []
    
    n = (end_interval - x0) / h 
    
    y = -1
    
    for i in range(int(n+1)):
        
        k_1 = h * function(x0,y0)
        k_2 = h * function(x0+(h/2) , y0 + (k_1/2))
        k_3 = h * function(x0 + (h/2), y0 + (k_2/2))
        k_4 = h * function((x0 + h), (y0 +k_3))
                    
        y = y0 + (k_1 / 6 + k_2 / 3 + k_3 / 3 + k_4 / 6)
        
        result.append((x0,y))
        
        y0 = y
        
        x0 = round(x0 + h, 5)

    print("The approximation of " + "y("  + str(end_interval) + ") with Runge-Kutta, " + "with h = " + str(h) + " is: " + str(result[int(n)][1]) + "\n")
    
    
    return result

data_e = runge_kutta(1,-1,0.05,2,f)

print("Assignment 1f): ")
print(linear_splines(data_e, Query_points, "Runge-Kutta"))

print("Exact solutions: \n")
print("1.052: " + str(exact_solution(1.052)))
print("1.555: " + str(exact_solution(1.555)))
print("1.978: " + str(exact_solution(1.978)) + "\n")

# Assignment 2
print("Assignment 2: \n")

def f2(x,y):
    return y - x ** 2

def exact_solution2(x):
    return 2 + 2 * x + x ** 2 - math.exp(x)


print("Exact solution y(3.3): " + str(exact_solution2(3.3)) + "\n")


def Adams_Bashforth(x0, y0, h, end_interval):
    
    '''
    Parameters
    ----------
    x0 : starting x value interval 
    y0 : starting y value
    h : step-size

    Returns
    -------
    result : interpolated values(x,y) in a list of tuples using Adams-Bashforth method.

    '''
    
    #Using exact solutions as starting values
    res = []
    res.append(exact_solution2(0))
    res.append(exact_solution2(0.1))
    res.append(exact_solution2(0.2))
    res.append(exact_solution2(0.3))
    
    n = end_interval / h
    result = []
    
    for i in range(len(res)):
        
        result.append((x0, res[i]))
        x0 = round(x0 + h, 5)
    
    # Now we start to iterate, x0 = 0.4 

    for i in range(3, int(n+1)):
        
        y = y0 + (h/24 * (55 * f2(x0-h,res[i]) - 59 * f2(x0-2*h, res[i-1]) + 37 * f2(x0 - 3 * h, res[i-2]) - 9 *f2(x0- 4 * h, res[i-3])))
        
        res.append(y)
        
        result.append((x0,y))
        
        y0 = y
        
        x0 = round(x0 + h, 5)
        
    
    print("The approximation of y(3.3) with Adams-Bashford, " + "with h = " + str(h) + " is: " + str(result[int(n+1)][1]) + "\n")
    
    
    return result


def Adams_Moulton(x0, y0, h, end_interval):
    '''

    Parameters
    ----------
    x0 : start interval x-coordinate
    y0 : starting y-value, which will be the third/fourth
    h : step-size
    end_interval : end interval x-coordinate

    Returns
    -------
    result : list of tuples (x,y) of interpolated values using Adams-Moulton.

    '''
    
    # Using exact solutions as starting values.
    res = []
    res.append(exact_solution2(0))
    res.append(exact_solution2(0.1))
    res.append(exact_solution2(0.2))
    res.append(exact_solution2(0.3))

    
    n = end_interval / h
    result = []
    
    # Adding initial values to the result list of tuples
    for i in range(len(res)):
        
        result.append((x0, res[i]))
        x0 = round(x0 + h, 5)
    
    # Start iteration, x0 = 0.4
    for i in range(3, int(n+1)):
        
        # Use Adams-Bashford first
        y_n = y0 + (h/24 * (55 * f2(x0-h,res[i]) - 59 * f2(x0-2*h, res[i-1]) + 37 * f2(x0 - 3 * h, res[i-2]) - 9 *f2(x0- 4 * h, res[i-3])))
        
        # Use the answer y_n above to correct value 
        y = y0 + (h/720 * (251 * f2(x0 ,y_n) + 646 * f2(x0-h, res[i]) - 264 * f2(x0 - 2 * h, res[i-1]) + 106 *f2(x0- 3 * h, res[i-2]) - 19 * f2(x0 - 4 * h, res[i-3])))
        
        res.append(y)
    
        result.append((x0,y))
        
        y0 = y
        
        x0 = round(x0 + h, 5)
        
    print("The approximation of y(3.3) with Adams-Moulton, " + "with h = " + str(h) + " is: " + str(result[int(n+1)][1]) + "\n")
    
    return result

Data_Runge = runge_kutta(0, 1, 0.1,3.3,f2)
data_AdamsBashforth = Adams_Bashforth(0, exact_solution2(0.3), 0.1,3.3)
data_AdamsMoulton = Adams_Moulton(0, exact_solution2(0.3), 0.1, 3.3)
