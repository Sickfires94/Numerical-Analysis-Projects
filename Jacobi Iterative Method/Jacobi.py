import numpy as np
import math
import random
from decimal import Decimal
import time

# Constants

DIAGONAL_VALUE: int = 1000
TOLERANCE_VALUE: int = 10 ** -20
TOLERANCE_CHECK: int = 5


# Return matrix of size n



def generateMatrix(size: int):
    matrix = []

    max_num = int(DIAGONAL_VALUE / (size - 1))
    for i in range(0, size):
        row = []

        for j in range(0, size):
            if i == j:
                row.insert(j, DIAGONAL_VALUE)
            else:
            	 row.insert(j, random.randint(0, max_num))

        matrix.insert(i, row)
    return matrix




# Returns the Result matrix of size n

def generateResult(size: int):
    result = []

    for i in range(0, size):
        result.insert(i, random.randint(1, 10))

    return result






# Returns Residual of Ax - B = 0

def Residual(matrix, coeff, result):
    residualMatrix = np.subtract(np.dot(matrix, coeff), result)
    sum = 0

    for i in range(0, len(coeff)):
        sum += Decimal(residualMatrix[i] ** 2)

    residual = math.sqrt(sum)

    return residual


# Performs Gauss Seidel method on A and B matrices and returns values of coefficients in an array

def Gauss_Seidel(matrix, result):
    coeff = []

    error = []

    for i in range(0, len(result)):
        coeff.insert(i, Decimal(0))

    counter = 0
    residual = Residual(matrix, coeff, result)
    error.append(residual)

    while residual > TOLERANCE_VALUE:

        for i in range(0, len(coeff)):

            sum: float = 0

            for j in range(0, len(coeff)):
                if (i == j):
                    continue
                sum += matrix[i][j] * coeff[j]

            coeff[i] = Decimal((result[i] - sum) / matrix[i][i])

        counter += 1

        if counter % TOLERANCE_CHECK == 0:
            residual = Residual(matrix, coeff, result)

    print("Gauss-Seidel needed iterations = ", counter)
    return coeff








# Performs Jacobi iterative method on A and B matrices and returns values of coefficients in an array

def Jacobi_iterative(matrix, result):
    coeff = []
    error = []

    for i in range(0, len(result)):
        value = Decimal(0)
        coeff.insert(i, value)

    counter = 0
    residual = Residual(matrix, coeff, result)
    error.insert(0, residual)

    while (residual > TOLERANCE_VALUE):
        currCoeff = coeff.copy()
        for i in range(0, len(currCoeff)):

            sum: float = 0

            for j in range(0, len(currCoeff)):
                if (i == j):
                    continue
                sum += matrix[i][j] * currCoeff[j]

            coeff[i] = Decimal((result[i] - sum) / matrix[i][i])

        counter += 1

        if counter % TOLERANCE_CHECK == 0:
            residual = Residual(matrix, coeff, result)


    print("Jacobi Iterative needed iterations = ", counter)
    return coeff



def main():
    size = 400
    matrix = generateMatrix(size)
    #print(matrix)

    result = generateResult(size)


    start_time = time.time()
    coeff = Gauss_Seidel(matrix, result)
    #print("Gauss Seidel co-efficients = ", coeff)

    print("time taken for Gauss Seidel Method = ","{:0.4f}".format(time.time() - start_time)," seconds")

    print("\n***************************************************** \n")

    start_time = time.time()
    coeff = Jacobi_iterative(matrix, result)
    #print("Jacobi iterative Co-efficients = ", coeff)
    print("time taken for Jacobi Iterative Method = ","{:0.4f}".format(time.time() - start_time)," seconds")
    print("\n***************************************************** \n")


main()
