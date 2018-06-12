
Program classical

Implicit None 
integer , parameter :: nx=351, ny=301,nn=10000,n=10000,m=10000
integer i,j, i1, j1,ii,jj
REAL*8 vh(nx,ny),vk(nx,ny),Ek(nx,ny),R(nx),P(ny)
REAL*8 x1(2),y1(2),z1(2),px1(2),py1(2),pz1(2),x2(2),y2(2),z2(2),px2(2),py2(2),pz2(2)&
&,U(2,2),k1(2),k2(2),k3(2),k4(2),EE1(n),EE2(n)&
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



!------------Ar--------------------
pai=3.1415926535897932384626433832795d0 
ee0=sqrt(1d15/3.51d16)
!w=0.057d0
W1=0.0584d0
W2=0.117d0
t0=2*pai/w1
!----------外场---------------------

n1=2.d0
n2=6.d0
h=0.1

!_______________________________________________________________________

do i=1,N

t1=h*i

enddo

!-------------------------------------------------------
open(file='xcos^2.dat',unit=15)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

xx1(1)=1.d0


pxx1(1)=0.d0



t(1)=0.d0



do i=1,1000

U(:,1)=(/xx1(1),pxx1(1)/)!!!CHUfHI


y11=xx1(1)
y12=pxx1(1)



t1=t(1)

K1=(/g1(y11,y12,t1),&
&f1(y11,y12,t1)/)


y21=xx1(1)+h*K1(1)/2
y22=pxx1(1)+h*K1(2)/2


t2=t(1)+h/2


K2=(/g1(y21,y22,t2),&
&f1(y21,y22,t2)/)

y31=xx1(1)+h*K2(1)/2
y32=Pxx1(1)+h*K2(2)/2



t3=t(1)+h/2

K3=(/g1(y31,y32,t3),&
&f1(y31,y32,t3)/)

y41=xx1(1)+h*K3(1)
y42=Pxx1(1)+h*K3(2)


t4=t(1)+h

K4=(/g1(y41,y42,t4),&
&f1(y41,y42,t4)/)

write(15,*)t1,xx1(1),dcos(t1)


do ii=1,2


U(ii,2)=U(ii,1)+(h*(K1(ii)+2.d0*K2(ii)+2.d0*K3(ii)+k4(ii))/6.d0)

end do



xx1(2)=U(1,2)
Pxx1(2)=U(2,2)




t(1)=t(1)+h


xx1(1)=xx1(2)
pxx1(1)=pxx1(2)


END DO                                         !时间循环


end


!------------------------函数---------------------------

FUNCTION f1(x1,px1,t1)   
Implicit None                                                 
REAL*8 f1,x1,Px1,t1

f1=-2 * sin(2 * t1) - 2 * x * cos(2 * t1)

return
END

 !!龙格库塔方法中用到的函数


FUNCTION  g1(x1,px1,t1) 
Implicit None 
REAL*8 x1,g1,px1,t1

g1=px1

return

end





!--------------------------------





