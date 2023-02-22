
import numpy as np

def neville(x_data, y_data, x):

  y = y_data.copy()
  m = len(x_data)  # number of data points

  for k in range(1, m):
    for i in range(m - k):
      y[i] = ((x - x_data[i+k])*y[i] + (x_data[i] - x)*y[i+1])/ \
      (x_data[i]-x_data[i+k])
      return (y[0])


# -----------------------------------------


def u_cal(u, n):

  temp = u
  for i in range(1, n):
    temp = temp * (u - i)
  return temp

def fact(n):
  f = 1
  for i in range(2, n + 1):
    f *= i
  return f
  
def newton_forward(x_data, y_data, approx):
  matrix = [[0 for i in range(len(x_data))] for j in range(len(x_data))]
  matrix[0][0] = y_data[0]
  matrix[1][0] = y_data[1]
  matrix[2][0] = y_data[2]
  matrix[3][0] = y_data[3]

  for i in range(1, len(x_data)):
    for j in range(len(x_data) - i):
      matrix[j][i] = matrix[j + 1][i - 1] - matrix[j][i - 1]
  sum = matrix[0][0]
  u = (approx - x_data[0]) / (x_data[1] - x_data[0])
  for i in range(1, len(x_data)):
    sum = sum + (u_cal(u, i) * matrix[0][i]) / fact(i)
  
  print(round(sum, 6))


#-------------------------------------------------
def hermite_approximation_matrix(x, fx, fdx):
    n = len(x)
    m = 2*n

    F = [[0]*m for i in range(n)]
    for i in range(n):
        F[i][0] = x[i]
        F[i][1] = fx[i]

    for j in range(2, m):
        for i in range(j - 1, n):
            if x[i] == x[i-j+1]:
                F[i][j] = fdx[i-j+1] / factorial(j-1)
            else:
                F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (x[i] - x[i-j+1])

    return F

#-------------------------------------------
def cubic_spline_interpolation(x, y):
  
  h = np.diff(x)
  u = np.zeros(len(x))
  v = np.zeros(len(x))
  b = np.zeros(len(x))
  n = len(x_data) - 1
  for i in range(1, len(x)-1):
      u[i] = 2*(h[i-1] + h[i])
      v[i] = 6*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

  # Create the matrix A
  A = np.zeros((len(x), len(x)))
  A[0, 0] = 1
  A[-1, -1] = 1
  for i in range(1, len(x)-1):
      A[i, i-1:i+2] = [h[i-1], u[i], h[i]]
      b[i] = 3 * (y[i+1] - y[i]) / h[i] - 3 * (y[i] - y[i-1]) / h[i-1]

  print(A)
  print("\n")
  print(b)
  u[0] = 0
  v[0] = b[0]

  for i in range(1, n):
    u[i] = h[i] / (2*(h[i-1]+h[i]) - h[i-1]*u[i-1])
    v[i] = (b[i] - h[i-1]*v[i-1]) / (2*(h[i-1]+h[i]) - h[i-1]*u[i-1])

  x = np.zeros(n+1)
  x[n] = v[n-1]

  for i in range(n-1, -1, -1):
    x[i] = v[i] - u[i]*x[i+1]

  print(x)


if __name__ == "__main__":

  # Question 1 : Neville's Method
  x_points = [3.6, 3.8, 3.9]
  y_points = [1.675, 1.436, 1.318]
  approximating_value = 3.7

  print(neville(x_points, y_points, approximating_value))

  #---------------------------------------

  # Question 2, 3 : Newton’s forward method 1, 2, 3 degrees

  x_data = [7.2, 7.4, 7.5, 7.6]
  y_data = [23.5492, 25.3913, 26.8224, 27.4589]

  appr_value = 7.3
  newton_forward(x_data, y_data, appr_value)

  #---------------------------------------
  # Question 4 :
  x = [3.6, 3.8, 3.9]
  fx = [1.675, 1.436, 1.318]
  fdx = [-1.195, -1.188, -1.182]
  F = hermite_approximation_matrix(x, fx, fdx)

  # print output
  for row in F:
      print(row)

  #---------------------------------------
  # Question 5 : cubic spline interpolation
  x1 = [2, 5, 8, 10]
  y1 = [3, 5, 7, 9]
  cubic_spline_interpolation(x1, y1)
