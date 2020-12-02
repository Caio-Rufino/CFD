import numpy as np
import Solver as S
import SolverImp as SI
import matplotlib.pyplot as pp
import time

nc = 150;
L = 10;
Ploc = 0.5;
dt = 0.01;
tf = 1;
v = 1.;
nt = int(tf/dt + 1);

tic = time.clock();
P = S.maccormack(L,v,Ploc,dt,nt,nc,6,6)
toc = time.clock();
time_b = toc - tic;

print(time_b)

lin, col = P.shape;

aux1 = np.zeros((lin));
aux1[:] = P[:,nt-1];

aux2 = np.zeros((lin));
aux2[:] = P[:,0];

l = np.linspace(0,L,lin);
t = np.linspace(0,tf,col);

titulo = 'MacCormack center6-center6';

pp.figure()
pp.contourf(t, l, P)
pp.title(titulo)
pp.xlabel('Tempo')
pp.ylabel('Local')
pp.savefig('IF_conv.png')

pp.figure()
pp.title(titulo)
pp.plot(l,aux2)
pp.plot(l,aux1)
pp.xlabel('Local')
pp.ylabel('Amplitude')
pp.legend(['Inicial','Final'])
pp.savefig('C_conv.png')

pp.show()
