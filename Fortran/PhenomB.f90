
!************************************************

!Title: Amplitude birefringence corrections in CS gravity

!Author: Kent Yagi

!Compile

! gfortran -c nrtype_qp.f90 (need to run only once)
! gfortran PhenomB.f90
! ./a.out

!Units: Everything is in second (in c=G=1 unit), e.g. 1M_sun = 5*10^-6 s
!************************************************


!Name of this program
program PhenomB

!we're going to use quadruple precision
use nrtype_qp

!This means that all the variables that you use need to be declared at first. Otherwise you'll get an error message.
implicit none


!Some of the variables may not be used in this program. don't worry too much.
integer, parameter :: imax=11, & !Number of parameters
nbin=100, & !number of bins for calculating the integral in the Fisher matrix
n1=2, & !1: GR, 2: non-GR
n_det=3, & !1: 1 detector (H), 2: 2 detectors (HL), 3: 3 detectors (HLV)
n_det_type=3, & !1: zero-detuned Adv. LIGO, 2: Adv. LIGO O1, 3:Cosmic Explorer (CE)
n_prior=1, & !1: with prior, 2: without prior
n_SNR_des=2 !1: desired SNR, 2: SNR determined from given luminosity distance
integer :: k,l,m,n
real(dp) :: m1, m2, mtotal, mc, eta, tsolar=4.925491_dp*10.0_dp**(-6), &
        dl, beta, tc, phic, tyear=3.1536_dp*10.0_dp**7, sigma, &
        fmax, fmin, tobs=1.0_dp, s=0.3_dp, &
        fend, fcut, fisco, fin,c=2.9979_dp*10.0_dp**10, &
        z,  d, snthr, integ, amp, kpc, dsub, sn,chi,f,f3,noiseligo, &
        lambdag,const,SNR_des,factor,sigma90,alpha,delta,psi,inc,GMST=0._dp, &
        GMST_H, GMST_L, GMST_V, Omega_s, phi_delay_H, phi_delay_L,phi_delay_V, &
        thetadot,H0,SNR,sigma_vel,n_gal,Omega_s_68,dc_min,dc_max,Vc_min,Vc_max,n_host,dl_rescale
real(dp), dimension(nbin) :: w, fbin
real(dp), dimension(imax,imax) :: fm, fmi,fm2,fmi2
complex(dpc) :: ad1,ad2,ad3,add1,add2

!************************************************
!Physical Parameters
!************************************************

!This is a factor that you need to multiply if you want 90% condfidence level instead of 1-sigma
sigma90=1.64485_dp


!Hubble constant
H0=70._dp*10.0_dp**5/(10.0_dp**6*3.261636_dp*tyear)/c

!Fiducial value for the time derivative of the scalar field
!Here, I'm assuming that GR is correct.
thetadot=0._dp

!************************************************
!Fiducial Values for Binary Parameters
!************************************************

    !desired SNR
    SNR_des=23.7_dp

    !redshift
    z=0.09_dp

    !luminosity distance
    dl=0.42_dp*10.0_dp**9*3.261636_dp*tyear

    !angle parameters
    alpha=8._dp*pi_d/24._dp
    delta=-75._dp*pi_d/180._dp
    psi=45._dp*pi_d/180._dp
    inc=60._dp*pi_d/180._dp

    !redshifted masses (in solar mass)
    m1=36.2_dp*(1.0_dp+z)
    m2=29.1_dp*(1.0_dp+z)

    !coalescence time and phase
    tc=0.0_dp
    phic=0.0_dp

    !effective spin
    chi=0.0_dp


!************************************************
!various masses
!************************************************

!total mass
!notice that I multiply this by tsolar to make it a unit in second
mtotal=(m1+m2)*tsolar

!symmetric mass ratio
eta=m1*m2/(mtotal/tsolar)**2

!chirp mass
mc=eta**(3.0_dp/5.0_dp)*mtotal



!************************************************
!waveform
!************************************************


!amplitude in the Fourier domain without the frequency dependence
const=(mtotal*eta**(3.0_dp/5.0_dp))**(5.0_dp/6.0_dp) &
/(pi_d**(2.0_dp/3.0_dp))/dl*sqrt(1.0_dp/30.0_dp)

!Calling the waveform function
call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,ad1,ad2,ad3,f3)

!************************************************
!Initial and final frequencies
!************************************************

!final frequency
fmax=f3

!Initial frequency
if (n_det_type==1 .or. n_det_type==2) then

    fmin=1.0_dp*10_dp**(1.0_dp)

else if (n_det_type==3) then

    fmin=1.0_dp

end if


!************************************************
!Fisher
!************************************************

call fish(mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,fm,fmi)

!90% credible in deg^2, cf. Eq.(13) in 1310.7454
Omega_s=-2._dp*pi_d*cos(delta)*log(1._dp-0.9_dp)*sqrt(fmi(7,7)*fmi(8,8)-fmi(7,8)*fmi(7,8))*3282.8_dp


print *, "S/N =",SNR
print *, "fmax=",f3
print *, "fmin=",fmin

print *, "delta mtotal/mtotal=", sqrt(fmi(1,1))
print *, "delta eta/eta=", sqrt(fmi(2,2))
print *, "delta chi=", sqrt(fmi(3,3))
print *, "delta dl/dl=", sqrt(fmi(4,4))
print *, "delta tc [s]=", sqrt(fmi(5,5))
print *, "delta phic [rad]=", sqrt(fmi(6,6))
print *, "delta alpha [rad]=", sqrt(fmi(7,7))
print *, "delta delta [rad]=", sqrt(fmi(8,8))
print *, "delta Omega_s [degree^2]=", Omega_s
print *, "delta psi [rad]=", sqrt(fmi(9,9))
print *, "delta inc [rad]=", sqrt(fmi(10,10))


print *, "delta DL [Mpc]=", sqrt(fmi(4,4))*dl/(10.0_dp**6*3.261636_dp*tyear)
!distance corresponding to the desired SNR
print *, "rescaled DL [Mpc]=", dl_rescale/(10.0_dp**6*3.261636_dp*tyear)


print *, "delta theta_dot=", sqrt(fmi(11,11))*c/(10._dp**5)

!measurement error in theta_dot when there is no correlation between parameters
print *, "delta theta_dot (uncor)=", 1._dp/sqrt(fm(11,11))*c/(10._dp**5)


!************************************************
!End of main part. Below are subroutines and functions
!************************************************



contains


!************************************************
!Subroutine for Fisher calculation
!************************************************

subroutine fish(mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,fm,fmi)

use nrtype_qp
implicit none

!input parameters
real(dp), intent(in) :: mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot

!output parameters (fisher matrix "fm" and its inverse "fmi")
real(dp), dimension(imax,imax), intent(inout) :: fm,fmi

integer :: i,j,n
real(dp) :: noise, amp, x, f, v, k4, a4, b4, c4, k5, a5, b5, c5, &
         snsa, s_gal, s_exgal, kappa=4.5_dp, dndf, shot, rad, accel, nn, &
         f0, clean=0.01_dp, noisedec,epsilon=10.0_dp**(-15), &
         aa0,aa1,aa2,aa3,aa4,aa5,aa6,lf,xf,s0,x1,p1,p2, &
         a1,a2,b1,b2,b3,b6,c1,c2,c3,noiseCE
complex(dpc) :: wi=(0.0_dp, 1.0_dp), adt1, adt2, adt3, ad1, ad2, adt_1, adt_2, adt_3, &
                add1, add2, addt1, addt2, addt_1, addt_2
complex(dpc), dimension(imax) ::dad1, dad2, dad3
real(dp), dimension(0,0) :: b
 

! frequency integration

!************************************************
!Initializing the Fisher matrix
!************************************************
!This is important when you're doing Monte Carlo simulation!
!It's always a good idea if you're adding things recursively. See my comments below.
fm=0.0_dp


!************************************************
!Calling Gauss-Legendre integrating subroutine
!************************************************

!This determins the frequency and weight for integration
call gauleg(fmin,fmax,fbin,w,nbin)



! do loop
do n=1, nbin

!renaming fbin to f
f=fbin(n)

!************************************************
!Taking derivative of the waveform with respect to each parameter using finite differencing
!************************************************


call wave(f,mtotal*(1._dp+epsilon),eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
!This second one is for complex conjugate
call wave(f,mtotal*(1._dp-epsilon),eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
!Hangford LIGO
dad1(1)=(adt1-adt_1)/(2.0_dp*epsilon)
!Livingston LIGO
dad2(1)=(adt2-adt_2)/(2.0_dp*epsilon)
!Virgo
dad3(1)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta*(1._dp+epsilon),chi,const,tc,phic,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta*(1._dp-epsilon),chi,const,tc,phic,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(2)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(2)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(2)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi+epsilon,const,tc,phic,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi-epsilon,const,tc,phic,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(3)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(3)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(3)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const*(1._dp+epsilon),tc,phic,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const*(1._dp-epsilon),tc,phic,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(4)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(4)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(4)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc+epsilon,phic,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc-epsilon,phic,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(5)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(5)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(5)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic+epsilon,alpha,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic-epsilon,alpha,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(6)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(6)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(6)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic,alpha+epsilon,delta,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic,alpha-epsilon,delta,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(7)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(7)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(7)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta+epsilon,psi,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta-epsilon,psi,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(8)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(8)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(8)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi+epsilon,inc,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi-epsilon,inc,thetadot,adt_1,adt_2,adt_3,f3)
dad1(9)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(9)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(9)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc+epsilon,thetadot,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc-epsilon,thetadot,adt_1,adt_2,adt_3,f3)
dad1(10)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(10)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(10)=(adt3-adt_3)/(2.0_dp*epsilon)

call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot+epsilon,adt1,adt2,adt3,f3)
call wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot-epsilon,adt_1,adt_2,adt_3,f3)
dad1(11)=(adt1-adt_1)/(2.0_dp*epsilon)
dad2(11)=(adt2-adt_2)/(2.0_dp*epsilon)
dad3(11)=(adt3-adt_3)/(2.0_dp*epsilon)


!************************************************
!Noise Curve
!************************************************


if (n_det_type==1) then

    !*******
    !LIGO zero-detuned by Ajith (2011)
    !*******

    f0=245.4_dp
    lf=f/f0

    !S_n(f)
    noise=10.0_dp**(-48.0_dp)*(0.0152*lf**(-4.0_dp)+0.2935*lf**(9.0_dp/4.0_dp)+2.7951_dp*lf**(3.0_dp/2.0_dp) &
    -6.508*lf**(3.0_dp/4.0_dp)+17.7622_dp)

else if (n_det_type==2) then

    !*******
    !LIGO O1
    !*******

    aa0=47.8466_dp
    aa1=-92.1896_dp
    aa2=35.9273_dp
    aa3=-7.61447_dp
    aa4=0.916742_dp
    aa5=-0.0588089_dp
    aa6=0.00156345_dp

    f0=log(f)

    !sqrt(S_n)
    noiseligo=sqrt(0.8464_dp)*exp(aa0 + aa1*f0 + aa2*f0**2 + aa3*f0**3 &
          + aa4*f0**4 + aa5*f0**5 + aa6*f0**6)

    !S_n(f)
    noise=noiseligo**2


else if (n_det_type==3) then

    !*******
    !CE
    !*******

    aa0=6.74257_dp
    aa1=-73.7462_dp
    aa2=34.509_dp
    aa3=-8.33831_dp
    aa4=1.09406_dp
    aa5=-0.073832_dp
    aa6=0.00201004_dp


    f0=log(f)

    noiseCE=exp(aa0 + aa1*f0 + aa2*f0**2 + aa3*f0**3 &
    + aa4*f0**4 + aa5*f0**5 + aa6*f0**6)

    noise=noiseCE**2


end if


!************************************************
!Fisher integration
!************************************************

!This is what I mean above by recursive summation.
!fm = fm + XXX

if (n_det_type==1 .or. n_det_type==2) then

    if (n_det==1) then

    do i=1,imax
    do j=1,imax
    fm(i,j)=fm(i,j) &
    +real((dad1(i)*conjg(dad1(j)))/noise) &
    *4.0_dp*w(n)
    end do
    end do

else if (n_det==2) then

    do i=1,imax
    do j=1,imax
    fm(i,j)=fm(i,j) &
    +real((dad1(i)*conjg(dad1(j))+dad2(i)*conjg(dad2(j)))/noise) &
    *4.0_dp*w(n)
    end do
    end do

else if (n_det==3) then

    do i=1,imax
    do j=1,imax
    fm(i,j)=fm(i,j) &
    +real((dad1(i)*conjg(dad1(j))+dad2(i)*conjg(dad2(j))+dad3(i)*conjg(dad3(j)))/noise) &
    *4.0_dp*w(n)
    end do
    end do

end if ! end n_det

else if (n_det_type==3) then

if (n_det==1) then

    do i=1,imax
    do j=1,imax
    fm(i,j)=fm(i,j) &
    +real((dad3(i)*conjg(dad3(j)))/noise) &
    *4.0_dp*w(n)
    end do
    end do

else if (n_det==3) then

    do i=1,imax
    do j=1,imax
    fm(i,j)=fm(i,j) &
    +real((dad1(i)*conjg(dad1(j))+dad2(i)*conjg(dad2(j))+dad3(i)*conjg(dad3(j)))/noise) &
    *4.0_dp*w(n)
    end do
    end do

end if ! end n_det

end if ! end n_det_type


!end of do loop
end do


!************************************************
!Rescaling SNR to the desired one
!************************************************


if (n_SNR_des==1) then

    dl_rescale=dl*sqrt(fm(4,4))/SNR_des

    fm=fm*SNR_des*SNR_des/fm(4,4)


else

    dl_rescale=dl

end if

!SNR
SNR=sqrt(fm(4,4))

!************************************************
!Prior on spin, coalescence phase and angles
!************************************************
if (n_prior==1 ) then

    fm(6,6)=fm(6,6)+(1.0_dp/(1._dp*pi_d))**2
    fm(3,3)=fm(3,3)+(1.0_dp/(1._dp))**2
    fm(7,7)=fm(7,7)+(1.0_dp/(pi_d/2._dp))**2
    fm(8,8)=fm(8,8)+(1.0_dp/(1._dp*pi_d))**2
    fm(9,9)=fm(9,9)+(1.0_dp/(pi_d/2._dp))**2
    fm(10,10)=fm(10,10)+(1.0_dp/(1._dp*pi_d))**2

end if


!************************************************
!Parameter estimation in GR
!(by setting the non-GR part of the Fisher matrix to a diagonal form)
!************************************************

if (n1==1) then

    do i=1,imax

    fm(i,11)=0.0_dp
    fm(11,i)=0.0_dp

    end do

    fm(11,11)=1.0_dp

end if


!************************************************
!Inverse of Fisher matrix
!************************************************


fmi=fm

!Calling a subroutine for calculating the inverse of a matrix
call gaussj(fmi,imax,imax,b,0,0)


end subroutine fish


!************************************************
!PhenomB waveform
!************************************************

subroutine wave(f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot,ad1,ad2,ad3,f3)

use nrtype_qp
implicit none

!input binary & theory parameters
real(dp), intent(in) :: f,mtotal,eta,chi,const,tc,phic,alpha,delta,psi,inc,thetadot

!output (final frequency and the waveform for H, L & V)
real(dp), intent(inout) :: f3
complex(dpc), intent(inout) :: ad1,ad2,ad3

integer :: i,j,n
real(dp) :: alpha2,alpha3,eps1,eps2,vf1,vf2,f1,f2,sigma,wm,wr,vf,phase,vCS
complex(dpc) :: wi=(0.0_dp, 1.0_dp)
real(dp), dimension(4,7) :: y
real(dp), dimension(5,7) :: z
real(dp), dimension(5) :: phi
real(dp), dimension(3) :: Fc,Fp,phi_pol

!************************************************
!Some parameters for PhenomB
!************************************************


alpha2=-323.0_dp/224.0_dp+451.0_dp*eta/168.0_dp
alpha3=(27.0_dp/8.0_dp-11.0_dp*eta/6.0_dp)*chi

eps1=1.4547_dp*chi-1.8897_dp
eps2=-1.8153_dp*chi+1.6557_dp


!************************************************
!Bunch of fitting parameters
!************************************************


y(1,1)=1.0_dp-4.455_dp*(1.0_dp-chi)**(0.217_dp)+3.521_dp*(1.0_dp-chi)**(0.26_dp)
y(2,1)=(1.0_dp-0.63_dp*(1.0_dp-chi)**(0.3_dp))/2.0_dp
y(3,1)=(1.0_dp-0.63_dp*(1.0_dp-chi)**(0.3_dp))*(1.0_dp-chi)**(0.45_dp)/4.0_dp
y(4,1)=0.3236_dp+0.04894_dp*chi+0.01346_dp*chi**(2.0_dp)

y(1,2)=0.6437_dp
y(1,3)=0.827_dp
y(1,4)=-0.2706_dp
y(1,5)=-0.05822_dp
y(1,6)=-3.935_dp
y(1,7)=-7.092_dp

y(2,2)=0.1469_dp
y(2,3)=-0.1228_dp
y(2,4)=-0.02609_dp
y(2,5)=-0.0249_dp
y(2,6)=0.1701_dp
y(2,7)=2.325_dp

y(3,2)=-0.4098_dp
y(3,3)=-0.03523_dp
y(3,4)=0.1008_dp
y(3,5)=1.829_dp
y(3,6)=-0.02017_dp
y(3,7)=-2.87_dp

y(4,2)=-0.1331_dp
y(4,3)=-0.08172_dp
y(4,4)=0.1451_dp
y(4,5)=-0.2714_dp
y(4,6)=0.1279_dp
y(4,7)=4.922_dp


z(1,1)=3715.0_dp/756.0_dp
z(2,1)=-16.0_dp*pi_d+113.0_dp*chi/3.0_dp
z(3,1)=15293365.0_dp/508032.0_dp-405.0_dp*chi**(2.0_dp)/8.0_dp
z(4,1)=0.0_dp
z(5,1)=0.0_dp

z(1,2)=-920.9_dp
z(1,3)=492.1_dp
z(1,4)=135.0_dp
z(1,5)=6742.0_dp
z(1,6)=-1053.0_dp
z(1,7)=-1.34_dp*10.0_dp**4

z(2,2)=1.702_dp*10.0_dp**4
z(2,3)=-9566.0_dp
z(2,4)=-2182.0_dp
z(2,5)=-1.214_dp*10.0_dp**5
z(2,6)=2.075_dp*10.0_dp**4
z(2,7)=2.386_dp*10.0_dp**5

z(3,2)=-1.254_dp*10.0_dp**5
z(3,3)=7.507_dp*10.0_dp**4
z(3,4)=1.338_dp*10.0_dp**4
z(3,5)=8.735_dp*10.0_dp**5
z(3,6)=-1.657_dp*10.0_dp**5
z(3,7)=-1.694_dp*10.0_dp**6

z(4,2)=-8.898_dp*10.0_dp**5
z(4,3)=6.31_dp*10.0_dp**5
z(4,4)=5.068_dp*10.0_dp**4
z(4,5)=5.981_dp*10.0_dp**6
z(4,6)=-1.415_dp*10.0_dp**6
z(4,7)=-1.128_dp*10.0_dp**7

z(5,2)=8.696_dp*10.0_dp**5
z(5,3)=-6.71_dp*10.0_dp**5
z(5,4)=-3.008_dp*10.0_dp**4
z(5,5)=-5.838_dp*10.0_dp**6
z(5,6)=1.514_dp*10.0_dp**6
z(5,7)=1.089_dp*10.0_dp**7

!************************************************
!Transition Frequencies
!************************************************

f1=(y(1,1)+y(1,2)*eta+y(1,3)*eta*chi+y(1,4)*eta*chi**2 &
    +y(1,5)*eta**2+y(1,6)*eta**2*chi+y(1,7)*eta**3)/(pi_d*mtotal)
f2=(y(2,1)+y(2,2)*eta+y(2,3)*eta*chi+y(2,4)*eta*chi**2 &
    +y(2,5)*eta**2+y(2,6)*eta**2*chi+y(2,7)*eta**3)/(pi_d*mtotal)
sigma=(y(3,1)+y(3,2)*eta+y(3,3)*eta*chi+y(3,4)*eta*chi**2 &
    +y(3,5)*eta**2+y(3,6)*eta**2*chi+y(3,7)*eta**3)/(pi_d*mtotal)
f3=(y(4,1)+y(4,2)*eta+y(4,3)*eta*chi+y(4,4)*eta*chi**2 &
    +y(4,5)*eta**2+y(4,6)*eta**2*chi+y(4,7)*eta**3)/(pi_d*mtotal)

!relative velocity at transition frequencies
vf1=(pi_d*mtotal*f1)**(1.0_dp/3.0_dp)
vf2=(pi_d*mtotal*f2)**(1.0_dp/3.0_dp)

!relative velocity at f
vf=(pi_d*mtotal*f)**(1.0_dp/3.0_dp)

!prefactor for merger and ringdown amplitude
wm=(1.0_dp+alpha2*vf1**2+alpha3*vf1**3)/(1.0_dp+eps1*vf1+eps2*vf1**2)
wr=pi_d*sigma/2.0_dp*(f2/f1)**(-2.0_dp/3.0_dp)*wm*(1.0_dp+eps1*vf2+eps2*vf2**2)

!************************************************
!Amplitude in GR
!************************************************

if (f < f1) then
        amp = const*f**(-7.0_dp/6.0_dp)*(1.0_dp+alpha2*vf**2+alpha3*vf**3)
        
  else if (f < f2) then
        amp = const*f1**(-7.0_dp/6.0_dp)*wm*(f/f1)**(-2.0_dp/3.0_dp) &
             *(1.0_dp+eps1*vf+eps2*vf**2)
       
  else if (f < f3) then 
        amp = const*f1**(-7.0_dp/6.0_dp)*wr*sigma/(2.0_dp*pi_d) &
             /((f-f2)**2+sigma**2/4.0_dp)
         
  else  
        amp = 0.0_dp
  
end if

!************************************************
!Phase in GR
!************************************************


do i=1,5

    phi(i)=z(i,1)+(z(i,2)+z(i,3)*chi+z(i,4)*chi**2)*eta+(z(i,5)+z(i,6)*chi)*eta**2+z(i,7)*eta**3

end do


phase=2.0_dp*pi_d*f*tc+phic+3.0_dp/(128.0_dp*eta*vf**5) &
    *(1.0_dp+phi(1)*vf**2+phi(2)*vf**3+phi(3)*vf**4+phi(4)*vf**6+phi(5)*vf**7)



!************************************************
!Beam-pattern functions for various detectors
!************************************************


!the arrival time at each detector
GMST_H=0.007204716433_dp*cos(delta)*cos(-1._dp*alpha+1._dp*GMST) &
-0.01278231726_dp*cos(delta)*sin(-1._dp*alpha+1._dp*GMST)-0.01533450076_dp*sin(delta)

GMST_L=0.0002475868047_dp*cos(delta)*cos(-1._dp*alpha+1._dp*GMST) &
-0.01832094573_dp*cos(delta)*sin(-1._dp*alpha+1._dp*GMST)-0.01074752339_dp*sin(delta)

GMST_V=-0.01515458033_dp*cos(delta)*cos(-1._dp*alpha+1._dp*GMST) &
+0.002809965657_dp*cos(delta)*sin(-1._dp*alpha+1._dp*GMST)-0.01459525654_dp*sin(delta)

!Arrival time delay among detectors
phi_delay_H=2._dp*pi_d*f*GMST_H
phi_delay_L=2._dp*pi_d*f*GMST_L
phi_delay_V=2._dp*pi_d*f*GMST_V


!Beam-pattern functions for Hanford (plus & cross modes)
Fp(1)=((-1.424276334_dp*cos(delta)**2+2.848552668_dp)*cos(psi)**2 &
-0.6209072974_dp*cos(psi)*sin(psi)*sin(delta)-1.424276335_dp &
+0.7121381674_dp*cos(delta)**2)*cos(-alpha+GMST_H)**2+(((-0.6209072976_dp &
+0.3104536479_dp*cos(delta)**2)*cos(psi)**2-2.848552669_dp*cos(psi)*sin(psi)*sin(delta) &
-0.1552268243_dp*cos(delta)**2+0.3104536479_dp)*sin(-alpha+GMST_H) &
+0.4947780891_dp*sin(delta)*cos(delta)-0.9119913319_dp*cos(psi)*sin(psi)*cos(delta) &
-0.9895561783_dp*cos(psi)**2*sin(delta)*cos(delta))*cos(-alpha+GMST_H) &
+(0.455995666_dp*sin(delta)*cos(delta)+0.9895561783_dp*cos(psi)*sin(psi)*cos(delta) &
-0.9119913318_dp*cos(psi)**2*sin(delta)*cos(delta))*sin(-alpha+GMST_H) &
+0.7121381671_dp+(-1.424276335_dp+0.4928680684_dp*cos(delta)**2)*cos(psi)**2 &
+0.3104536479_dp*cos(psi)*sin(psi)*sin(delta)-0.2464340343_dp*cos(delta)**2


Fc(1)=(-0.6209072975_dp*cos(psi)**2*sin(delta)+(-2.848552669_dp*sin(psi) &
+1.424276336_dp*sin(psi)*cos(delta)**2)*cos(psi)+0.3104536479_dp*sin(delta))*cos(-alpha+GMST_H)**2 &
+((-2.848552668_dp*cos(psi)**2*sin(delta)+(0.620907297_dp*sin(psi) &
-0.310453648_dp*sin(psi)*cos(delta)**2)*cos(psi)+1.424276334_dp*sin(delta))*sin(-alpha+GMST_H) &
-0.9119913318_dp*cos(psi)**2*cos(delta)+0.9895561784_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta) &
+0.4559956657_dp*cos(delta))*cos(-alpha+GMST_H)+(-0.4947780891_dp*cos(delta) &
+0.9895561781_dp*cos(psi)**2*cos(delta) &
+0.9119913318_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta))*sin(-alpha+GMST_H) &
+0.310453648_dp*cos(psi)**2*sin(delta)+(1.424276335_dp*sin(psi) &
-0.4928680684_dp*sin(psi)*cos(delta)**2)*cos(psi)-0.1552268242_dp*sin(delta)

!Beam-pattern functions for Livingston
Fp(2)=((1.040573092_dp*cos(delta)**2-2.081146184_dp)*cos(psi)**2 &
+1.121682131_dp*cos(psi)*sin(psi)*sin(delta)+1.040573093_dp &
-0.5202865458_dp*cos(delta)**2)*cos(-alpha+GMST_L)**2 &
+(((1.12168213_dp-0.5608410653_dp*cos(delta)**2)*cos(psi)**2 &
+2.081146184_dp*cos(psi)*sin(psi)*sin(delta) &
+0.2804205326_dp*cos(delta)**2-0.5608410653_dp)*sin(-alpha+GMST_L) &
-0.4945891839_dp*sin(delta)*cos(delta)+0.72646252_dp*cos(psi)*sin(psi)*cos(delta) &
+0.9891783677_dp*cos(psi)**2*sin(delta)*cos(delta))*cos(-alpha+GMST_L) &
+(-0.36323126_dp*sin(delta)*cos(delta)-0.9891783677_dp*cos(psi)*sin(psi)*cos(delta) &
+0.7264625205_dp*cos(psi)**2*sin(delta)*cos(delta))*sin(-alpha+GMST_L) &
-0.5202865457_dp+(1.040573093_dp+0.3865389484_dp*cos(delta)**2)*cos(psi)**2 &
-0.5608410652_dp*cos(psi)*sin(psi)*sin(delta)-0.193269474_dp*cos(delta)**2

Fc(2)=(1.12168213_dp*cos(psi)**2*sin(delta)+(2.081146185_dp*sin(psi) &
-1.040573091_dp*sin(psi)*cos(delta)**2)*cos(psi)-0.5608410652_dp*sin(delta))*cos(-alpha+GMST_L)**2 &
+((2.081146184_dp*cos(psi)**2*sin(delta)+(-1.121682131_dp*sin(psi) &
+0.5608410653_dp*sin(psi)*cos(delta)**2)*cos(psi) &
-1.040573093_dp*sin(delta))*sin(-alpha+GMST_L)+0.7264625204_dp*cos(psi)**2*cos(delta) &
-0.9891783677_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta) &
-0.3632312599_dp*cos(delta))*cos(-alpha+GMST_L)+(0.494589184_dp*cos(delta) &
-0.9891783677_dp*cos(psi)**2*cos(delta) &
-0.7264625204_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta))*sin(-alpha+GMST_L) &
-0.5608410652_dp*cos(psi)**2*sin(delta)+(-1.040573093_dp*sin(psi) &
-0.3865389482_dp*sin(psi)*cos(delta)**2)*cos(psi)+0.2804205326_dp*sin(delta)


!Beam-pattern functions for Virgo
Fp(3)=((1.383399753_dp*cos(delta)**2-2.766799504_dp)*cos(psi)**2 &
-0.7926702315_dp*cos(psi)*sin(psi)*sin(delta)+1.383399752_dp &
-0.6916998761_dp*cos(delta)**2)*cos(-alpha+GMST_V)**2+(((-0.7926702315_dp &
+0.3963351156_dp*cos(delta)**2)*cos(psi)**2+2.766799504_dp*cos(psi)*sin(psi)*sin(delta) &
-0.198167558_dp*cos(delta)**2+0.3963351156_dp)*sin(-alpha+GMST_V) &
+0.465152434_dp*sin(delta)*cos(delta)-0.7513324117_dp*cos(psi)*sin(psi)*cos(delta) &
-0.9303048681_dp*cos(psi)**2*sin(delta)*cos(delta))*cos(-alpha+GMST_V) &
+(0.3756662059_dp*sin(delta)*cos(delta)+0.9303048681_dp*cos(psi)*sin(psi)*cos(delta) &
-0.7513324117_dp*cos(psi)**2*sin(delta)*cos(delta))*sin(-alpha+GMST_V) &
-0.6916998761_dp+(1.383399752_dp-1.303555288_dp*cos(delta)**2)*cos(psi)**2 &
+0.3963351156_dp*cos(psi)*sin(psi)*sin(delta)+0.6517776447_dp*cos(delta)**2


Fc(3)=(-0.7926702315_dp*cos(psi)**2*sin(delta) &
+(2.766799505_dp*sin(psi)-1.383399754_dp*sin(psi)*cos(delta)**2)*cos(psi) &
+0.3963351156_dp*sin(delta))*cos(-alpha+GMST_V)**2+((2.766799504_dp*cos(psi)**2*sin(delta) &
+(0.7926702316_dp*sin(psi)-0.3963351156_dp*sin(psi)*cos(delta)**2)*cos(psi) &
-1.383399752_dp*sin(delta))*sin(-alpha+GMST_V)-0.7513324117_dp*cos(psi)**2*cos(delta) &
+0.9303048681_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta) &
+0.3756662059_dp*cos(delta))*cos(-alpha+GMST_V)+(-0.465152434_dp*cos(delta) &
+0.9303048681_dp*cos(psi)**2*cos(delta) &
+0.7513324117_dp*cos(psi)*sin(psi)*sin(delta)*cos(delta))*sin(-alpha+GMST_V) &
+0.3963351156_dp*cos(psi)**2*sin(delta)+(-1.383399753_dp*sin(psi) &
+1.303555289_dp*sin(psi)*cos(delta)**2)*cos(psi)-0.198167558_dp*sin(delta)



!parity violation parameter
vCS=-pi_d*f*H0*dl*thetadot


!************************************************
!Polarization phase
!************************************************


phi_pol(1)=atan(2._dp*cos(inc)*Fc(1)/((1._dp+cos(inc)*cos(inc))*Fp(1))) &
+vCS*((1._dp-cos(inc)*cos(inc))**2*Fp(1)*Fc(1)) &
/((1._dp+cos(inc)*cos(inc))**2*Fp(1)*Fp(1)+4._dp*cos(inc)*cos(inc)*Fc(1)*Fc(1))
phi_pol(2)=atan(2._dp*cos(inc)*Fc(2)/((1._dp+cos(inc)*cos(inc))*Fp(2))) &
+vCS*((1._dp-cos(inc)*cos(inc))**2*Fp(2)*Fc(2)) &
/((1._dp+cos(inc)*cos(inc))**2*Fp(2)*Fp(2)+4._dp*cos(inc)*cos(inc)*Fc(2)*Fc(2))
phi_pol(3)=atan(2._dp*cos(inc)*Fc(3)/((1._dp+cos(inc)*cos(inc))*Fp(3))) &
+vCS*((1._dp-cos(inc)*cos(inc))**2*Fp(3)*Fc(3)) &
/((1._dp+cos(inc)*cos(inc))**2*Fp(3)*Fp(3)+4._dp*cos(inc)*cos(inc)*Fc(3)*Fc(3))


!************************************************
!Final Waveform (with amplitude correction included)
!************************************************

ad1=amp*5._dp/4._dp*(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(1)*Fp(1)+4._dp*cos(inc)*cos(inc)*Fc(1)*Fc(1)) &
    +vCS*(2._dp*cos(inc)*(1._dp+cos(inc)*cos(inc))*(Fp(1)*Fp(1)+Fc(1)*Fc(1))) &
    /(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(1)*Fp(1)+4._dp*cos(inc)*cos(inc)*Fc(1)*Fc(1)))) &
    *exp(wi*(phase-phi_pol(1)+phi_delay_H))
ad2=amp*5._dp/4._dp*(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(2)*Fp(2)+4._dp*cos(inc)*cos(inc)*Fc(2)*Fc(2)) &
    +vCS*(2._dp*cos(inc)*(1._dp+cos(inc)*cos(inc))*(Fp(2)*Fp(2)+Fc(2)*Fc(2))) &
    /(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(2)*Fp(2)+4._dp*cos(inc)*cos(inc)*Fc(2)*Fc(2)))) &
    *exp(wi*(phase-phi_pol(2)+phi_delay_L))
ad3=amp*5._dp/4._dp*(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(3)*Fp(3)+4._dp*cos(inc)*cos(inc)*Fc(3)*Fc(3)) &
    +vCS*(2._dp*cos(inc)*(1._dp+cos(inc)*cos(inc))*(Fp(3)*Fp(3)+Fc(3)*Fc(3))) &
    /(sqrt((1._dp+cos(inc)*cos(inc))**2*Fp(3)*Fp(3)+4._dp*cos(inc)*cos(inc)*Fc(3)*Fc(3)))) &
    *exp(wi*(phase-phi_pol(3)+phi_delay_V))



end subroutine wave



!************************************************
!Gauss-Legendre integration subroutine
!(taken from Numerical Recipe)
!************************************************

SUBROUTINE gauleg(x1,x2,x,w,n)

INTEGER :: i,j,m,n
REAL(DP) :: x1,x2,EPS=3.q-14,p1,p2,p3,pp,xl,xm,z,z1
REAL(DP), DIMENSION(n) :: x,w

m=(n+1)/2
xm=0.5q0*(x2+x1)
xl=0.5q0*(x2-x1)
do 12 i=1,m
z=cos(3.141592654q0*(i-.25q0)/(n+.5q0))
1       continue
p1=1.q0
p2=0.q0
do 11 j=1,n
p3=p2
p2=p1
p1=((2.q0*j-1.q0)*z*p2-(j-1.q0)*p3)/j
11        continue
pp=n*(z*p1-p2)/(z*z-1.q0)
z1=z
z=z1-p1/pp
if(abs(z-z1).gt.EPS)goto 1
x(i)=xm-xl*z
x(n+1-i)=xm+xl*z
w(i)=2.q0*xl/((1.q0-z*z)*pp*pp)
w(n+1-i)=w(i)
12    continue
return
END SUBROUTINE



!************************************************
!inverse matrix subroutine
!(taken from Numerical Recipe)
!************************************************

SUBROUTINE gaussj(a,n,np,b,m,mp)

use nrtype_qp
implicit none

INTEGER :: m,mp,n,np
INTEGER, PARAMETER :: NMAX=50
REAL(DP), DIMENSION(np,np) :: a
REAL(DP), DIMENSION(np,mp) :: b
INTEGER :: i,icol,irow,j,k,l,ll
INTEGER, DIMENSION(NMAX) :: indxc,indxr,ipiv
REAL(DP) :: big,dum,pivinv
do 11 j=1,n
ipiv(j)=0
11    continue
do 22 i=1,n
big=0.q0
do 13 j=1,n
if(ipiv(j).ne.1)then
do 12 k=1,n
if (ipiv(k).eq.0) then
if (abs(a(j,k)).ge.big)then
big=abs(a(j,k))
irow=j
icol=k
endif
else if (ipiv(k).gt.1) then
pause 'singular matrix in gaussj'
endif
12          continue
endif
13      continue
ipiv(icol)=ipiv(icol)+1
if (irow.ne.icol) then
do 14 l=1,n
dum=a(irow,l)
a(irow,l)=a(icol,l)
a(icol,l)=dum
14        continue
do 15 l=1,m
dum=b(irow,l)
b(irow,l)=b(icol,l)
b(icol,l)=dum
15        continue
endif
indxr(i)=irow
indxc(i)=icol
if (a(icol,icol).eq.0.q0) pause 'singular matrix in gaussj'
pivinv=1.q0/a(icol,icol)
a(icol,icol)=1.q0
do 16 l=1,n
a(icol,l)=a(icol,l)*pivinv
16      continue
do 17 l=1,m
b(icol,l)=b(icol,l)*pivinv
17      continue
do 21 ll=1,n
if(ll.ne.icol)then
dum=a(ll,icol)
a(ll,icol)=0.q0
do 18 l=1,n
a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
do 19 l=1,m
b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
endif
21      continue
22    continue
do 24 l=n,1,-1
if(indxr(l).ne.indxc(l))then
do 23 k=1,n
dum=a(k,indxr(l))
a(k,indxr(l))=a(k,indxc(l))
a(k,indxc(l))=dum
23        continue
endif
24    continue
return

END subroutine gaussj

end program PhenomB
