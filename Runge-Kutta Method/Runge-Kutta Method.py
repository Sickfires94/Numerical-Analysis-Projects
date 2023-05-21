import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import math

HEIGHT = 1  # Original Set Height
LOWERLIMIT = 0 
UPPERLIMIT = 87 # a = last two digits of my ERP = 87 
INITIAL_Y = 87/1.23

constants = [1,2,3,4]
theta = [0.25,0.5,0.75]
weights = [0.25,0.5,0.75]

INITIAL_X = LOWERLIMIT

def ComputeDerivitive(y,x):
    dydt = (87*x + 1.23*y)/(1.23*x + 87*y)
    return dydt

def CorrectAnswer(x):
    t = [0,x]
    y = odeint(ComputeDerivitive,INITIAL_Y,t)
    return y[1]

def computeAbsoluteSquareError(x,y):
    correct_answer = CorrectAnswer(x)
    return (y - correct_answer) ** 2

def constant_sum(constants):
    sum = 0
    for i in constants:
        sum += i
    return sum
    
def RK_Method(h, initial_y, initial_x):
    error = 0
    y = initial_y
    x = initial_x
    constant =  constant_sum(constants)
    for i in range(int(UPPERLIMIT/h)):
        
        k1 = h * ComputeDerivitive(y,x)
        k2 = h * ComputeDerivitive(y + k1 * theta[0],x + h * weights[0])
        k3 = h * ComputeDerivitive(y + k2 * theta[1], x + h * weights[1])
        k4 = h * ComputeDerivitive(y + k3 * theta[2], x + h * weights[2])
        
        y = y + ((constants[0] * k1 + constants[1] * k2 + constants[2] * k3 + constants[3] * k4) / constant)
        x += h
        #print(computeAbsoluteSquareError(x,y))
        error += computeAbsoluteSquareError(x,y)
     
    return math.sqrt((error)/(UPPERLIMIT/h))
        
def main():
    errors = []
    heights = []
    for i in range(0,5):
        errors.append(RK_Method(HEIGHT/(2 ** i), INITIAL_Y, INITIAL_X))
        print(" RMSE = ", errors[-1], ", at Height = ", HEIGHT / 2 ** i)
        heights.append(HEIGHT/(2 ** i))

    print(heights)

    plt.plot(heights, errors)
    plt.xlabel("Height")
    plt.ylabel("RMSE")
    plt.show()
        
main()
