import pandas as pd
import numpy as np
import scipy
import control
from statsmodels.tsa.vector_ar.var_model import VAR
import statsmodels.api as sm

class OptimalControl:

    def __init__(self):
        pass
 
    def getB(self, yaw, deltat):
        B = np.array([[np.cos(yaw)*deltat, 0],[np.sin(yaw)*deltat, 0],[0, deltat]])
        return B

    def state_space_model(self, A, state_t_minus_1, B, control_input_t_minus_1):
        state_estimate_t = (A @ state_t_minus_1) + (B @ control_input_t_minus_1) 
        return state_estimate_t
    
    def lqr(self, actual_state_x, desired_state_xf, Q, R, A, B, dt):
        x_error = actual_state_x - desired_state_xf
        N = 1
        P = [None] * (N + 1)
        Qf = Q
        P[N] = Qf
        for i in range(N, 0, -1):
            P[i-1] = Q + A.T @ P[i] @ A - (A.T @ P[i] @ B) @ np.linalg.pinv(R + B.T @ P[i] @ B) @ (B.T @ P[i] @ A)      
        K = [None] * N
        u = [None] * N
        for i in range(N):
            K[i] = -np.linalg.pinv(R + B.T @ P[i+1] @ B) @ B.T @ P[i+1] @ A
            u[i] = K[i] @ x_error
            u_star = u[N-1]
        return u_star