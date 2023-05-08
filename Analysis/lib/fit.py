import numpy as np

def gauss(x,amp,mean,std):
    return amp*np.exp(-np.power(x - mean, 2.) / (2 * np.power(std, 2.)))

def scint_profile(x,const,a_f,tau_f,tau_s):
    return const*(2*a_f/tau_f*np.exp(-(x)/tau_f) + 2*(1-a_f)/tau_s*np.exp(-(x)/tau_s))