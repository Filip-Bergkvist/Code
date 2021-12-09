import matplotlib.pyplot as plt
import math

datan = [(1,2),(2,2.5),(3,7),(4,10.5),(6,12.75),(8,13),(10,13)]

# X coordinates 
x_query = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5]


def linear_splines(data,x_q):
    '''

    Parameters
    ----------
    data : list of tuples with (x,y) values 
    x_q : query points (x coordinates)

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
    
    plt.title("Linear")
    plt.plot(x_val, y_val)
    plt.show()
    return "Done \n"

print(linear_splines(datan, x_query))


def interpolation(data, value):
    '''
    

    Parameters
    ----------
    data : list of tupled with(x,y) values
    value : query point

    Returns
    -------
    sm : interpolated value for the query point

    '''
    y_value = 0
    
    for k in range(len(data)):
            
        t = 1
        
        for j in range(len(data)):
            
            if k != j:
                
                
                t = t *  ((value - data[j][0]) / (data[k][0] - data[j][0]))
            
        
        y_value = y_value + (data[k][1] * t)
        
    
    return y_value



def lagrange_interpolate(data, x_q):
    '''
    
    Parameters
    ----------
    data : list of tuples with (x,y) values
    x_q : query points (x coordinates)

    Returns
    -------
    Interpolated values for the query points using lagrange.

    '''
    result = []

    for i in x_q:
        
        result.append((i,interpolation(data, i)))

    print("Interpolated values Lagrange: \n")
    for i in range(len(result)):
       print("X coordinate: " + " " + str(result[i][0]) + " " + "Interpolated value: " + str(result[i][1])+ "\n")
        
    x_val = [x[0] for x in result]
    y_val = [x[1] for x in result]
    
    plt.title("Lagrange")
    plt.plot(x_val, y_val)
    plt.show()
        
    return "Done"


print(lagrange_interpolate(datan, x_query))

#Funktion from lab assignment:
def f(x,y):
    
    return y - x

print("Euler method: \n")

def euler_method(x0,y0,h):
    
    result = []
    
    n = 1 / h
    
    y_prim = 0
    
    for i in range(int(n+1)):
        
        y_prim = f(x0,y0)
        
        y = y0 + h * y_prim
        
        result.append((x0,y))
        
        x0 = round(x0 + h, 5)
    
        y0 = y

    
    print("The approximation of y(1) " + "with h = " + str(h) + " is: " + str(result[int(n)][1]) + "\n")
    
    return result[int(n)][1]
    
    
# Present results

figure = []

exact_solution = 1 + 1 - math.exp(1) / 2

figure.append(("h = 0.1",(euler_method(0,0.5,0.1))))
figure.append(("h = 0.05",(euler_method(0,0.5,0.05))))
figure.append(("h = 0.01",(euler_method(0,0.5,0.01))))
figure.append(("Exact solution", exact_solution))

for q,w in figure: 
    
    print(f"{q}\t\t{w}")

