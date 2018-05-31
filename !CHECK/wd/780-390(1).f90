
Program classical

Implicit None 
integer , parameter :: nx=351, ny=301,nn=10000,n=10000,m=10000
integer i,j, i1, j1,ii,jj
REAL*8 vh(nx,ny),vk(nx,ny),Ek(nx,ny),R(nx),P(ny)
REAL*8 x1(2),y1(2),z1(2),px1(2),py1(2),pz1(2),x2(2),y2(2),z2(2),px2(2),py2(2),pz2(2)&
&,U(12,2),k1(12),k2(12),k3(12),k4(12),EE1(n),EE2(n)&
&,rx1(m),ry1(m),rz1(m),rpx1(m),rpy1(m),rpz1(m),rx2(m),ry2(m),rz2(m),rpx2(m),rpy2(m),rpz2(m)&
&,xx1(2),yy1(2),zz1(2),pxx1(2),pyy1(2),pzz1(2),xx2(2),yy2(2),zz2(2),pxx2(2),pyy2(2),pzz2(2)&
&,xxa(n),yya(n),zza(n),pxxa(n),pyya(n),pzza(n),xxb(n),yyb(n),zzb(n),pxxb(n),pyyb(n),pzzb(n),rrr(n)

REAL*8 a,q,E0,Pai,RR,PP,aa2,aa3,aa4,aa5,t1,t2,t3,t4,t(1),aa,bb,tion,&
&y11,y12,y13,y14,y15,y16,y17,y18,y19,y111,y112,y113,y114,&
&y21,y22,y23,y24,y25,y26,y27,y28,y29,y211,y212,y213,y214,&
&y31,y32,y33,y34,y35,y36,y37,y38,y39,y311,y312,y313,y314,&
&y41,y42,y43,y44,y45,y46,y47,y48,y49,y411,y412,y413,y414,E(n)
REAL*8 ee0,n1,n2,h,w,t0,fwhm,qq1,tdi,tre,trec,aaa1,za,zb,tao,tp,tt,c,lamda,W1,W2,zz
Real*8 matrix(nx, ny), x(nx), y(ny), min,z,pd,ip2,s,pion
Real*8 p1, p2, r1, r2,r12,ef
REAL*8,EXTERNAL::g1,g2,g3,g4,g5,g6,f1,f2,f3,f4,f5,f6,E1,E2,QQ



open(file='interaction.DAT',unit=20)
open(file='E.DAT',unit=21)
open(file='Ei-12.DAT',unit=22)
open(file='pd.DAT',unit=23)
open(file='sdi.DAT',unit=24)
open(file='nsdi.DAT',unit=25)
open(file='ob4.DAT',unit=30)
open(file='test-n.DAT',unit=300)
!-------------------------------
open(file='k1.DAT',unit=200)
open(file='k2.DAT',unit=201)
open(file='k3.DAT',unit=202)
open(file='k4.DAT',unit=203)
open(file='ob.DAT',unit=208)
open(file='wd_k1_k2_k3_k4.dat',unit=210)
!------------Ar--------------------
 pai=3.1415926535897932384626433832795d0 
 a=2.d0
 q=1.225d0
 E0=-1.59d0
 Ip2=-1.065d0
!----------�ⳡ---------------------

ee0=sqrt(1d15/3.51d16)
!w=0.057d0
W1=0.0584d0
W2=0.117d0
t0=2*pai/w1
n1=2.d0
n2=6.d0
!h=(2.d0*n1+n2)*t0/dble(n)
h=1

tp=0.8d0
!-------------------------------
!c=2.997925d0*(10.d0**8.d0)
!lamda=740.d0*(10.d0**(-9.d0)) 
!w1=2.d0*pai*c/lamda         
!w1=w1*2.4189d0*(10.d0**(-17.d0))
!t0=2.d0*pai/w1
!fwhm=7.d0*(10.d0**(-15.d0))                         !!! ��߿� �� (����ǿ��)��
!tao=fwhm/2.d0/( sqrt(2.d0*log(sqrt(2.d0)))  )       !!! tao �룻 f(t)=E0 * exp[ -1/2 (t/tao)^2 ]
!tao=tao/(2.4189d0*10.d0**(-17.d0)) 
!tt=15.d0*(10.d0**(-15.d0))                          !!! ��
!Tp=tt/(2.4189d0*(10.d0**(-17.d0)))                  !!! ԭ�ӵ�λ
!h=2.d0*Tp/(n-1)
!n2=0.7d0                                           !��ƫ��
!EE0=2.742*10**3.d0*sqrt(4d15)    
!EE0=EE0/(5.1421*(10**11.D0))
!n1=0.d0*pai
!------------------------------------
do i=1,N

t1=h*i
E(i)=sqrt(E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0+E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp)**2.d0)
write(21,'(1X,21(1X,F15.7))')t1/t0,E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp),E(i)

enddo

!_______________________________________________________________________

open(file='Initialization.dat',unit=15)

do j=1,m

read(15,*)rx1(j),ry1(j),rz1(j),rpx1(j),rpy1(j),rpz1(j),rx2(j),ry2(j),rz2(j),rpx2(j),rpy2(j),rpz2(j)

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do jj=21,21,1

EE0=2.742*10**3.d0*sqrt((10.d0**(12.d0+dble(jj-1)*0.2d0)) )
!EE0=2.742*10**3.d0*sqrt((jj*2.d0)*(10**13.D0))    
EE0=EE0/(5.1421*(10**11.D0))
ee0=0.d0

!--------------------------------
zz=0.d0
z=0.d0



do j=1,1



CALL RANDOM_NUMBER(aa5)

!tao=2.d0*aa5*pai

xx1(1)=rx1(j)
xx2(1)=rx2(j)
yy1(1)=ry1(j)
yy2(1)=ry2(j)
zz1(1)=rz1(j)
zz2(1)=rz2(j)

pxx1(1)=rpx1(j)
pxx2(1)=rpx2(j)
pyy1(1)=rpy1(j)
pyy2(1)=rpy2(j)
pzz1(1)=rpz1(j)
pzz2(1)=rpz2(j)


aa=0.d0
bb=0.d0


t(1)=0.d0



do i=1,n

U(:,1)=(/xx1(1),Pxx1(1),xx2(1),Pxx2(1),yy1(1),pyy1(1),yy2(1),pyy2(1),zz1(1),pzz1(1),zz2(1),pzz2(1)/)!!!CHUfHI


y11=xx1(1)
y12=Pxx1(1)
y13=xx2(1)
y14=Pxx2(1)
y15=yy1(1)
y16=Pyy1(1)
y17=yy2(1)
y18=Pyy2(1)
y111=zz1(1)
y112=pzz1(1)
y113=zz2(1)
y114=pzz2(1)


t1=t(1)

K1=(/g1(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f1(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g2(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f2(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&g3(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f3(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q)-E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g4(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f4(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q)-E2(t1,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g5(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f5(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q)-E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g6(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q),&
&f6(y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114,t1,a,q)-E1(t1,t0,ee0,w1,w2,n1,n2,tao,tp)/)


y21=xx1(1)+h*K1(1)/2
y22=pxx1(1)+h*K1(2)/2
y23=xx2(1)+h*K1(3)/2
y24=pxx2(1)+h*K1(4)/2
y25=yy1(1)+h*K1(5)/2
y26=Pyy1(1)+h*K1(6)/2
y27=yy2(1)+h*K1(7)/2
y28=Pyy2(1)+h*K1(8)/2
y211=zz1(1)+h*K1(9)/2
y212=pzz1(1)+h*K1(10)/2
y213=zz2(1)+h*K1(11)/2
y214=pzz2(1)+h*K1(12)/2

t2=t(1)+h/2


K2=(/g1(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f1(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g2(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f2(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&g3(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f3(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q)-E2(t2,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g4(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f4(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q)-E2(t2,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g5(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f5(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q)-E1(t2,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g6(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q),&
&f6(y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214,t2,a,q)-E1(t2,t0,ee0,w1,w2,n1,n2,tao,tp)/)

y31=xx1(1)+h*K2(1)/2
y32=Pxx1(1)+h*K2(2)/2
y33=xx2(1)+h*K2(3)/2
y34=Pxx2(1)+h*K2(4)/2
y35=yy1(1)+h*K2(5)/2
y36=Pyy1(1)+h*K2(6)/2
y37=yy2(1)+h*K2(7)/2
y38=Pyy2(1)+h*K2(8)/2
y311=zz1(1)+h*K2(9)/2
y312=Pzz1(1)+h*K2(10)/2
y313=zz2(1)+h*K2(11)/2
y314=Pzz2(1)+h*K2(12)/2



t3=t(1)+h/2

K3=(/g1(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f1(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g2(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f2(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&g3(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f3(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q)-E2(t3,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g4(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f4(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q)-E2(t3,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g5(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f5(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q)-E1(t3,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g6(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q),&
&f6(y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314,t3,a,q)-E1(t3,t0,ee0,w1,w2,n1,n2,tao,tp)/)

y41=xx1(1)+h*K3(1)
y42=Pxx1(1)+h*K3(2)
y43=xx2(1)+h*K3(3)
y44=Pxx2(1)+h*K3(4)
y45=yy1(1)+h*K3(5)
y46=Pyy1(1)+h*K3(6)
y47=yy2(1)+h*K3(7)
y48=Pyy2(1)+h*K3(8)
y411=zz1(1)+h*K3(9)
y412=Pzz1(1)+h*K3(10)
y413=zz2(1)+h*K3(11)
y414=Pzz2(1)+h*K3(12)

t4=t(1)+h

K4=(/g1(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f1(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g2(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f2(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&g3(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f3(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q)-E2(t4,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g4(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f4(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q)-E2(t4,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g5(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f5(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q)-E1(t4,t0,ee0,w1,w2,n1,n2,tao,tp),&
&g6(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q),&
&f6(y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414,t4,a,q)-E1(t4,t0,ee0,w1,w2,n1,n2,tao,tp)/)

do ii=1,12
U(ii,2)=U(ii,1)+(h*(K1(ii)+2.d0*K2(ii)+2.d0*K3(ii)+k4(ii))/6.d0)

end do

write(210,'(1X,21(1X,F15.12))')y11,y21,y31,y41
write(210,'(1X,21(1X,F15.12))')y12,y22,y32,y42
write(210,'(1X,21(1X,F15.12))')y13,y23,y33,y43
write(210,'(1X,21(1X,F15.12))')y14,y24,y34,y44
write(210,'(1X,21(1X,F15.12))')y15,y25,y35,y45
write(210,'(1X,21(1X,F15.12))')y16,y26,y36,y46
write(210,'(1X,21(1X,F15.12))')y17,y27,y37,y47
write(210,'(1X,21(1X,F15.12))')y18,y28,y38,y48
write(210,'(1X,21(1X,F15.12))')y111,y211,y311,y411
write(210,'(1X,21(1X,F15.12))')y112,y212,y312,y412
write(210,'(1X,21(1X,F15.12))')y113,y213,y313,y413
write(210,'(1X,21(1X,F15.12))')y114,y214,y314,y414



!write(200,'(1X,21(1X,F15.12))')y11,y12,y13,y14,y15,y16,y17,y18,y111,y112,y113,y114
!write(201,'(1X,21(1X,F15.12))')y21,y22,y23,y24,y25,y26,y27,y28,y211,y212,y213,y214
!write(202,'(1X,21(1X,F15.12))')y31,y32,y33,y34,y35,y36,y37,y38,y311,y312,y313,y314
!write(203,'(1X,21(1X,F15.12))')y41,y42,y43,y44,y45,y46,y47,y48,y411,y412,y413,y414









xx1(2)=U(1,2)
Pxx1(2)=U(2,2)
xx2(2)=U(3,2)
Pxx2(2)=U(4,2)
yy1(2)=U(5,2)
Pyy1(2)=U(6,2)
yy2(2)=U(7,2)
Pyy2(2)=U(8,2)
zz1(2)=U(9,2)
Pzz1(2)=U(10,2)
zz2(2)=U(11,2)
Pzz2(2)=U(12,2)

!write(208,'(1X,21(1X,F15.12))')t(1),xx1(2),yy1(2),zz1(2),pxx1(2),pyy1(2),pzz1(2),xx2(2),yy2(2),zz2(2),pxx2(2),pyy2(2),pzz2(2)

 
EE1(i)=pxx1(2)**2.d0/2.d0+pyy1(2)**2.d0/2.d0+pzz1(2)**2.d0/2.d0-2.d0/sqrt(xx1(2)**2.d0+yy1(2)**2.d0+zz1(2)**2.d0)&
   &+1.d0/(xx1(2)**2.d0+yy1(2)**2.d0+zz1(2)**2.d0)*q**2.d0/4.d0/a*exp(a*(1.d0-((xx1(2)**2.d0+yy1(2)**2.d0+zz1(2)**2.d0)&
   &*(pxx1(2)**2.d0+pyy1(2)**2.d0+pzz1(2)**2.d0)/q**2.d0)**2.d0))&
   &+1.d0/sqrt((xx1(2)-xx2(2))**2.d0+(yy1(2)-yy2(2))**2.d0+(zz1(2)-zz2(2))**2.d0)/2.d0
              

EE2(i)=pxx2(2)**2.d0/2.d0+pyy2(2)**2.d0/2.d0+pzz2(2)**2.d0/2.d0-2.d0/sqrt(xx2(2)**2.d0+yy2(2)**2.d0+zz2(2)**2.d0)&
   &+1.d0/(xx2(2)**2.d0+yy2(2)**2.d0+zz2(2)**2.d0)*q**2.d0/4.d0/a*exp(a*(1.d0-((xx2(2)**2.d0+yy2(2)**2.d0+zz2(2)**2.d0)&
   &*(pxx2(2)**2.d0+pyy2(2)**2.d0+pzz2(2)**2.d0)/q**2.d0)**2.d0))&
   &+1.d0/sqrt((xx1(2)-xx2(2))**2.d0+(yy1(2)-yy2(2))**2.d0+(zz1(2)-zz2(2))**2.d0)/2.d0

rrr(i)=sqrt((xx1(2)-xx2(2))**2.d0+(yy1(2)-yy2(2))**2.d0+(zz1(2)-zz2(2))**2.d0)

t(1)=t(1)+h

xx1(1)=xx1(2)
pxx1(1)=pxx1(2)
xx2(1)=xx2(2)
pxx2(1)=pxx2(2)          !!���¸�ֵ   ����
yy1(1)=yy1(2)
pyy1(1)=pyy1(2)
yy2(1)=yy2(2)
pyy2(1)=pyy2(2)
zz1(1)=zz1(2)
pzz1(1)=pzz1(2)
zz2(1)=zz2(2)
pzz2(1)=pzz2(2)



xxa(i)=xx1(2)
pxxa(i)=pxx1(2)
xxb(i)=xx2(2)
pxxb(i)=pxx2(2)          !!���¸�ֵ   ����
yya(i)=yy1(2)
pyya(i)=pyy1(2)
yyb(i)=yy2(2)
pyyb(i)=pyy2(2)
zza(i)=zz1(2)
pzza(i)=pzz1(2)
zzb(i)=zz2(2)
pzzb(i)=pzz2(2)




END DO                  !ʱ��ѭ��
!write(300,'(1X,21(1X,F15.12))')xx1(2),yy1(2),zz1(2),pxx1(2),pyy1(2),pzz1(2),xx2(2),yy2(2),zz2(2),pxx2(2),pyy2(2),pzz2(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(ee1(n)*ee2(n).lt.0.d0)then

zz=zz+1.d0

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if((ee1(n).gt.0.d0).and.(ee2(n).gt.0.d0))then

write(20,'(1X,21(1X,F15.7))')real(j),pxx1(2),pxx2(2),pyy1(2),pyy2(2),pzz1(2),pzz2(2)


z=z+1.d0




!!!!!!!!!!!!!!!!SDI-NSDI!!!!!!!!!!!!!!!!!!!

!---------������ʱ��--------------------

do ii=1,n

if(EE1(ii)*EE2(ii).lt.0.d0)then

tion=ii

exit
end if
end do


!--------˫����ʱ��--------
do ii=1,n

if((EE1(ii).gt.0.d0).and.(EE2(ii).gt.0.d0))then

tdi=ii

exit
end if
end do




!-------------------------
aaa1=0.d0

do i=tion,tdi
if(rrr(i).gt.5.d0)then
tre=i
exit
end if
end do


do i=tre,tdi
if(rrr(i).lt.5.d0)then
aaa1=1.d0
trec=i
exit
end if
end do




if(aaa1.eq.1.d0)then

write(25,'(1X,21(1X,F15.7))')real(j),pxx1(2),pxx2(2),pyy1(2),pyy2(2),pzz1(2),pzz2(2),tion*h/t0,trec*h/t0,tdi*h/t0

zb=zb+1.d0

end if

if(aaa1.eq.0.d0)then

write(24,'(1X,21(1X,F15.7))')real(j),pxx1(2),pxx2(2),pyy1(2),pyy2(2),pzz1(2),pzz2(2),tion*h/t0,tdi*h/t0

za=za+1.d0
end if

!!!!!!!!!!!!!һ����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if(EE1(tion).gt.0.d0)then
   
aa=EE1(tion)
bb=e0-aa
else
aa=EE2(tion)
bb=e0-aa
write(22,*)bb,aa


end if


end if

end do    !������ѭ��



pd=z/m
pion=zz/m
s=zb/za

!-----------------------------------

write(23,'(1X,21(1X,F15.7))')real(jj),pd,pion,s

end do



end

!-----------------------------------------------------------
subroutine find_max(nx, ny, matrix, min, i, j)
Implicit None 
integer i, j,i1,j1,nx,ny
real*8 matrix(nx, ny), min,tem 


min=matrix(1,1)
i=1
j=1
do j1=1, ny
	do i1=1, nx
	if(min.ge.matrix(i1,j1)) then
	min=matrix(i1,j1)
	i=i1
	j=j1
	
	
	end if
	end do
end do

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------����---------------------------

FUNCTION f1(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q)   
Implicit None                                                  !!��������������õ��ĺ���
REAL*8 f1,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f1=x1*  ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
   &+(x1-x2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)




return
END


FUNCTION f2(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None                  
REAL*8 f2,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q




f2=x2*  ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(x1-x2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)



return
END

FUNCTION f3(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q)   
Implicit None 
REAL*8 f3,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f3=y1* ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
   &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
    &+(y1-y2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)




return
END



FUNCTION f4(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 f4,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q

f4=y2* ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
    &*exp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(y1-y2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)


return
end


FUNCTION f5(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None  
REAL*8 f5,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f5=z1* ( (q**2.d0/2.d0/a/(z1**2.d0+x1**2.d0+y1**2.d0)**2.d0+(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0)&
    &*exp(a*(1.d0-((z1**2.d0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0) )&
    &+(z1-z2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)



return
END



FUNCTION f6(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 f6,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q


f6=z2* ( (q**2.d0/2.d0/a/(z2**2.d0+x2**2.d0+y2**2.d0)**2.d0+(pz2**2.d0+px2**2.d0+py2**2.d0)**2.d0/q**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0))-2.d0/sqrt((z2**2.d0+y2**2.d0+x2**2.d0)**3.d0) )&
    &-(z1-z2)/sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
!---------------


return
end




FUNCTION  g1(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g1,a,q

g1=px1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
     &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )
  

return

end


FUNCTION  g2(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g2,a,q

g2=px2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )


return

end



FUNCTION  g3(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g3,a,q

g3=py1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
      &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )

return
end


FUNCTION  g4(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g4,a,q

g4=py2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
      &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )

return

end


FUNCTION  g5(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g5,a,q

g5=pz1*( 1.d0 - 1.d0/q**2.D0*(z1**2.D0+y1**2.d0+x1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)&
    &*dexp(a*(1.d0-((z1**2.d0+x1**2.d0+y1**2.d0)*(pz1**2.d0+px1**2.d0+py1**2.d0)/q**2.d0)**2.d0)) )

return

end


FUNCTION  g6(x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,t,a,q) 
Implicit None 
REAL*8 t,x1,Px1,x2,Px2,y1,Py1,y2,Py2,z1,Pz1,z2,Pz2,g6,a,q

g6=pz2*( 1.d0 - 1.d0/q**2.D0*(z2**2.D0+y2**2.d0+x2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)&
    &*dexp(a*(1.d0-((z2**2.d0+x2**2.d0+y2**2.d0)*(pz2**2.d0+px2**2.d0+py2**2.d0)/q**2.d0)**2.d0)) )

return
end 

!-------------�ⳡ����------------------------------------------------------



Function  E1(t,t0,ee0,w1,w2,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E1,QQ,h,fwhm,n,tao,tp,w1,w2


!-------------T------------------------                                              !ԭ�ӵ�λ
!if(t.lt.n1*T0) then
!QQ=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!QQ=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!QQ=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!QQ=0.D0
!endif
!----------sin2-------------------------
QQ=(DSIN(w1/2.D0/(2*n1+n2)*t))**2.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!E1=EE0*QQ*(sin(w1*t)-sin(w2*t))


E1=(ee0/(1.d0+tp))*QQ*(sin(w1*t+tao))-(ee0*tp/(1.d0+tp))*QQ*(sin(w2*t+2*tao))
!-------------ellipse---------------


!E1=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 )*(n2/sqrt(1.d0+n2**2.d0))*SIN(w1*t+N1)
!-------------------------------------



!E1=ee0*QQ*sin(w1*t)
return
end

Function  E2(t,t0,ee0,w1,w2,n1,n2,tao,tp)
Implicit None 
REAL*8 ee0,w,t,n1,n2,t0,E2,QQ,h,fwhm,n,tao,tp,w1,w2
  

!-------------T------------------------                                            
!if(t.lt.n1*T0) then
!QQ=t/(n1*T0)
!end if
!if((t.ge.n1*T0).and.(t.le.(n1+n2)*T0)) then
!QQ=1.d0
!endif
!if((t.gt.(n1+n2)*T0).and.(t.le.(2.d0*n1+n2)*T0)) then
!QQ=-t/(n1*T0)+(2.d0*n1+n2)/n1
!endif
!if(t.GT.(2.d0*n1+n2)*T0) then
!QQ=0.D0
!endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!E2=EE0*QQ*(COS(w1*t)+COS(w2*t))

!----------sin2-------------------------
QQ=(DSIN(w1/2.D0/(2*n1+n2)*t))**2.D0


E2=(ee0/(1.d0+tp))*QQ*(cos(w1*t+tao))+(ee0*tp/(1.d0+tp))*QQ*(cos(w2*t+2*tao))

!-------------ellipse---------------

!E2=ee0*exp( -((t-tp)/tao)**2.d0/2.d0 ) *(1.d0/sqrt(1.d0+n2**2.d0))*cos(w1*t+N1)

!------------------------------------------

!E2=ee0*QQ*cos(w1*t)
return
end

!--------------------------------





