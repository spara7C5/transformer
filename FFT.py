####################################################
######  PSD: direct vs library methods comparison ##
####################################################

from pylab import *
from scipy import signal
from scipy import random

from FUNDEF_rand import RandomWalk
########-----entries-------######



filtering=0
mirroring=0

######## opt 1) reading data #################
#tvec,Svec=genfromtxt("../data",delimiter='\t',usecols=(0,1),max_rows=1000000,unpack=True)
#Ts=(tvec[11]-tvec[10])
#lo=cos(2*pi*(10**6)*tvec)
#S=Svec*lo
##############################################

myrand=RandomWalk()
######## opt 2) generating noise #############
Ts=1*10**(-3) 
Fs=1/Ts
psd_noise= 7 # [V^2/Hz]
amp=sqrt(psd_noise*(Fs/2)) # [V]
######### 2_1) white noise##############
#Svec=random.normal(0,sqrt(7*500),10000)
########################################

########## 2_2 rand walk ##############
amprand=sqrt((10**6)*2*pi*pi/Fs)
Svec=[(amprand)*myrand.funrand(1) for i in range(100000) ]
#######################################

tvec=array([Ts*i for i in range(100000)])
#fig0=figure();plot(tvec,Svec);title("rand noise")
S=Svec
#######################################

Fs=1/Ts
N=len(Svec)


if mirroring:
	S=pad(S,(0,len(S)),'symmetric')



#########filtering################

if filtering :
	fc=70 #cut frequency
	b,a=signal.butter(3, fc/(Fs/2), 'low')
	Sf=signal.lfilter(b, a, S)
	Sfl=signal.filtfilt(b, a, S) #phaselinearfilt

	fig1=figure()
	w, h = signal.freqz(b, a)
	plot(w*(Fs/(2*pi)), abs(h))
	grid()
	Ysf = fft(Sf)
	Ysfl = fft(Sfl)
	Ysfmod,Ysflmod =abs(Ysf),abs(Ysfl)


####----FFT---------########

fig2=figure()
ax2=fig2.add_subplot(111)
Ys = fft(S,len(S)) #length is half of data-length
Ysmod=2*((abs(Ys)**2)/(Fs*N**1))
fvec=(arange(len(S)))*(Fs/len(S)) 
ax2.plot(fvec[0:int(N/2)],(Ysmod)[0:int(N/2)])
#ax2.ticklabel_format(axis='x',style='sci',scilimits=(0,6) )
xax=ax2.get_xaxis().get_major_formatter()
xax.set_powerlimits((1,6))
xax.set_scientific(True)
ax2.grid(linestyle="--")
#####IMPORTANT: f=0 must be omitted in PSD ##########
print(mean(Ysmod[1:int(N/2)]))
ax2.semilogy()


fpsd,Ppsd=signal.welch(S,fs=Fs,window='boxcar',scaling='density',nperseg=len(S))
ax2.plot(fpsd,Ppsd,c=(0,1,0))
print(mean(Ppsd))
#fper,Pper=signal.periodogram(S,fs=Fs)
#plot(fper,Pper,c=(1,0,0))
#ax3=fig2.add_subplot(212)
#psd(S,5000,Fs=Fs)

show()



