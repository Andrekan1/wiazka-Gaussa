# -*- coding: utf-8 -*-
# Uruchomienie: python wykresy_dyspersja.py   (wymaga numpy i matplotlib)
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
# Odtwarza wykresy z prezentacja/img/ na podstawie danych wyeksportowanych z COMSOL (ten folder).
HERE=os.path.dirname(os.path.abspath(__file__))
SRC=HERE+'/'
OUT=os.path.join(HERE,'..','img')+'/'
SURF='#fcfcfb'; T1='#0b0b0b'; T2='#52514e'; GRID='#e4e3df'; REF='#8a8984'
C=['#2a78d6','#eb6834','#1baf7a','#eda100']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.edgecolor':GRID,'axes.labelcolor':T2,
  'xtick.color':T2,'ytick.color':T2,'axes.facecolor':SURF,'figure.facecolor':SURF,'axes.spines.top':False,'axes.spines.right':False})
lam=1e-6; w0=10e-6; h=3e-6; k0=2*np.pi/lam
def eta(n,z): return np.arctan(z/(np.pi*n*w0**2/lam))
def wrap(a): return (a+np.pi)%(2*np.pi)-np.pi
d=np.genfromtxt(SRC+'dyspersja_numeryczna.csv',delimiter=',',comments='%')
R={}
for n,kf,Nl,_,ph,amp in d: R.setdefault((n,kf),[]).append((Nl,ph,amp))
def get(n,kf):
    r=sorted(R[(n,kf)]); N=np.array([a[0] for a in r]); ph=np.array([a[1] for a in r]); A=np.array([a[2] for a in r])
    return N,np.abs(wrap(ph-eta(n,h))),A
def style(ax):
    ax.grid(True,which='major',color=GRID,lw=0.8); ax.set_axisbelow(True)
def head(fig,t,s):
    fig.suptitle(t,x=0.06,y=0.975,ha='left',fontsize=21,fontweight='bold',color=T1)
    fig.text(0.06,0.915,s,ha='left',fontsize=13,color=T2)
# ---------- Fig 1: phase error
cases=[((1.0,1.0),'k₁ = n·k₀ (dopasowany): n = 1',C[0],'o','-'),
       ((1.5,1.5),'k₁ = n·k₀ (dopasowany): n = 1,5',C[0],'o',':'),
       ((1.5,1.0),'k₁ ≠ n·k₀: n = 1,5, k₁ = k₀',C[1],'s','-'),
       ((1.0,1.5),'k₁ ≠ n·k₀: n = 1, k₁ = 1,5·k₀',C[2],'D','-'),
       ((1.0,0.0),'k₁ = 0 (pełne pole): n = 1',C[3],'^','-')]
fig,ax=plt.subplots(figsize=(16,9),dpi=100)
for key,lab,col,mk,ls in cases:
    N,err,_=get(*key)
    mfc=col if ls=='-' else SURF
    ax.plot(N,err,color=col,lw=2.5,ls=ls,marker=mk,ms=11,mfc=mfc,mec=col if ls!='-' else SURF,mew=2,label=lab,zorder=3)
ax.annotate('k₁ = n·k₀ (n = 1 i n = 1,5 – krzywe się pokrywają)',(8,get(1.0,1.0)[1][-1]),xytext=(12,0),textcoords='offset points',va='center',color=T1,fontsize=13)
ax.annotate('k₁ ≠ n·k₀ (n = 1,5; k₁ = k₀)',(8,get(1.5,1.0)[1][-1]),xytext=(12,-8),textcoords='offset points',va='center',color=T1,fontsize=13)
ax.annotate('k₁ ≠ n·k₀ (n = 1; k₁ = 1,5·k₀)',(8,get(1.0,1.5)[1][-1]),xytext=(12,8),textcoords='offset points',va='center',color=T1,fontsize=13)
ax.annotate('k₁ = 0 – pełne pole (n = 1)',(8,get(1.0,0.0)[1][-1]),xytext=(12,0),textcoords='offset points',va='center',color=T1,fontsize=13)
ax.set_xscale('log',base=2); ax.set_yscale('log')
ax.set_xticks([1,2,4,8]); ax.set_xticklabels(['1','2','4','8'])
ax.set_xlim(0.85,40); ax.set_ylim(3e-4,8)
ax.axhline(np.pi,color=REF,ls='--',lw=1.2); ax.text(0.9,np.pi*1.15,'π – faza całkowicie przypadkowa',color=T2,fontsize=12)
style(ax)
ax.set_xlabel('Nl – liczba elementów wzdłuż z na h = 3 µm (= 3λ)')
ax.set_ylabel('|błąd fazy| na osi przy z = h  [rad]')
ax.legend(loc='lower left',bbox_to_anchor=(0.01,0.12),frameon=False,fontsize=12)
head(fig,'Dyspersja numeryczna: błąd fazy wiązki HG(0,0) po 3 µm propagacji',
     'błąd = arg(Ex·exp(+i·n·k₀·h)) − faza Gouya z teorii;  COMSOL ewbe, siatka poprzeczna „Normal”, sweep n, k₁, Nl')
fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(OUT+'dyspersja_blad_fazy.png'); plt.close(fig)
# ---------- Fig 2: amplitude
fig,ax=plt.subplots(figsize=(16,9),dpi=100)
amp_cases=[((1.0,1.0),'k₁ = n·k₀: n = 1',C[0],'o','-'),((1.5,1.5),'k₁ = n·k₀: n = 1,5',C[0],'o',':'),
           ((1.5,1.0),'n = 1,5, k₁ = k₀',C[1],'s','-'),((1.0,1.5),'n = 1, k₁ = 1,5·k₀',C[2],'D','-')]
for key,lab,col,mk,ls in amp_cases:
    N,_,A=get(*key); mfc=col if ls=='-' else SURF
    ax.plot(N,A,color=col,lw=2.5,ls=ls,marker=mk,ms=11,mfc=mfc,mec=col if ls!='-' else SURF,mew=2,label=lab,zorder=3)
for val,txt in [(1.0,'1 – teoria (ośrodek bezstratny)'),(0.8,'2k₁/(k₁+n·k₀) = 0,8'),(1.2,'2k₁/(k₁+n·k₀) = 1,2')]:
    ax.axhline(val,color=REF,ls='--',lw=1.2,zorder=1); ax.text(9.0,val+0.012,txt,va='bottom',color=T2,fontsize=13)
ax.set_xscale('log',base=2); ax.set_xticks([1,2,4,8]); ax.set_xticklabels(['1','2','4','8']); ax.set_xlim(0.85,40); ax.set_ylim(0.2,1.4)
style(ax); ax.set_xlabel('Nl – liczba elementów wzdłuż z'); ax.set_ylabel('|Ex| / E₀ na osi przy z = h')
ax.legend(loc='lower right',bbox_to_anchor=(0.62,0.02),frameon=False,fontsize=12)
head(fig,'Amplituda na wyjściu: niedopasowane k₁ psuje też warunek brzegowy',
     'przy k₁ ≠ n·k₀ amplituda dąży do 2k₁/(k₁+n·k₀) zamiast do 1 – zagęszczanie siatki tego nie naprawia')
fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(OUT+'dyspersja_amplituda.png'); plt.close(fig)
# ---------- field along axis
f=np.genfromtxt(SRC+'pole_na_osi.txt',comments='%')
def series(n,kf,Nl):
    m=(f[:,1]==n)&(f[:,2]==kf)&(f[:,3]==Nl); return f[m,0],f[m,5]
def figfield(n,sel,fname,title,sub):
    fig,ax=plt.subplots(figsize=(16,9),dpi=100)
    for (kf,Nl,lab,col,ls) in sel:
        zz,re=series(n,kf,Nl); ax.plot(zz*1e6,re,color=col,lw=2.5,ls=ls,label=lab,zorder=3)
    z=np.linspace(0,h,601); wz=w0*np.sqrt(1+(z/(np.pi*n*w0**2/lam))**2)
    ax.plot(z*1e6,(w0/wz)*np.cos(n*k0*z-eta(n,z)),color=T1,lw=1.6,ls=(0,(4,3)),label='teoria (wiązka Gaussa)',zorder=4)
    ax.set_xlabel('z [µm]'); ax.set_ylabel('Re(Ex) na osi [V/m]'); ax.set_xlim(0,3); ax.set_ylim(-1.4,1.4)
    style(ax); ax.legend(loc='lower left',frameon=False,fontsize=13,ncol=2)
    head(fig,title,sub); fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(OUT+fname); plt.close(fig)
figfield(1.0,[(1,1,'obwiednia k₁ = n·k₀ = k₀, Nl = 1',C[0],'-'),(0,4,'pełne pole k₁ = 0, Nl = 4 (1,3 el./λ)',C[3],':'),(0,8,'pełne pole k₁ = 0, Nl = 8 (2,7 el./λ)',C[3],'-')],
  'dyspersja_pole_n1.png','n = 1: obwiednia trafia w teorię już przy 1 elemencie, pełne pole – nie',
  'pełne pole przy < 3 elementach na λ: zła długość fali (dyspersja numeryczna) i zła amplituda')
figfield(1.5,[(1.5,1,'k₁ = n·k₀, Nl = 1',C[0],'-'),(1,2,'k₁ = k₀, Nl = 2',C[1],':'),(1,8,'k₁ = k₀, Nl = 8',C[1],'-')],
  'dyspersja_pole_n15.png','n = 1,5: wybór k₁ decyduje, ile elementów trzeba',
  'k₁ = k₀ → obwiednia niesie exp(−i·0,5·k₀z) i potrzebuje gęstej siatki;  k₁ = n·k₀ → 1 element wystarcza')
print('ok')
