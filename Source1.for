!C**********************************************************************
!C               SUBROUTINE UMAT ¨C SuperDislocation (SD) Model          *
!C                    Version: Wagoner 1.0.2 (08/10/2018)              *
!C                  Case: Tension test of a sample with                *
!C                      102,040  elements and 15 grains                *
!C======================================================================
!C| The SD Model code was developed by the Wagoner Group at the        * 
!C| Department of Materials Science and Engineering, The Ohio State    *
!C| University.  Developed starting in 2004 from a crystal             *
!C| plasticity (CP) program of Kalidindi, SD UMAT versions were created* 
!C| sequentially by  Myoung-Gyu Lee, Hojun Lim, Hyuk Jong Bong,        *
!C| Jisheng Qin & R. H. Wagoner.                                       *
!C|                                                                    *
!C| Version 1.01 was developed primarily from Lim 5.7.3 and secondarily* 
!C| from Hyuk1 and Hyuk2 versions.  The major changes by Jisheng Qin   *
!C| are as follows:                                                    *
!C| O Module ¡°modname¡± with dynamic memory allocation replaced most of  *
!C|   the previous Common Block / static memory storage.               *
!C| O Output of the NYE/GND tensor was corrected to reflect distinct   *
!C|   K13 and K31 components.                                          *
!C| O ¡°NSTORE,¡± an integer array, was previously and erroneously        *
!C|   declared as Real*8 in UMAT.  This was corrected.                 *
!C| O The REDIST and REDIST1 subroutines, for Lim and Hyuk versions,   *
!C|   respectively, were replaced by REDIST3. REDIST3 computes the     *
!C|   global displacement vector R from target element (local, stress) * 
!C|   to a remote element (Superdislocation).  This vector is expressed*
!C|   in local coordinates to get X1, X2, X3, where X1 is slip/burgers *
!C|   vector direction, X2 is the slip plane normal, and X3 is in the  *
!C|   X1 x X2 (line) direction.                                        *
!C**********************************************************************
     
      Module modname

cccccccccccccccccccccccccccccccc
c-twq
      CHARACTER(79) FILEroutine
	Include 'para.txt'

c      PARAMETER(MAXNOEL=196,MAXNPT=8,MAXSLIP=24)
c	 PARAMETER(MAXGRAIN=4,MAXROW=2000,Maxcandi=1000)
c      PARAMETER(maxnode=100000,maxkount=100)
!     ASSIGNS # TO UNITS AND OPENS I/O FILES.
C	
      PARAMETER(UR0=15)     ! SD.IN 
      PARAMETER(UR1=16)     ! ABAQUS.INP  
      PARAMETER(UR2=17)     ! abaqus_info.dat 	
	PARAMETER(UR3=101)    ! tex_testa.OUT 
      PARAMETER(UR4=102)    ! coord_testa.OUT
      PARAMETER(UR5=103)    ! astrain.OUT
	PARAMETER(UR6=104)    ! lstrain.OUT
	PARAMETER(UR7=105)    ! meg.OUT
      PARAMETER(UR8=106)    ! check_testa.OUT  
      PARAMETER(UR9=107)    ! iter.OUT
      PARAMETER(UR10=108)   ! DATA.OUT 
	PARAMETER(UR11=109)   ! transmissivity.out 
     
      parameter ( PI=4.*ATAN(1.))

	CHARACTER(25) FILEINP

      INTEGER INIT,INCR
      real(8)  ONE(3,3),ZERO(3,3),
     1      NCRYS,NSLIP,SMATC(MAXSLIP,3,3),
     2   	  FPITALL(MAXNOEL,MAXNPT,3,3),
     4   	  TSTARTALL(MAXNOEL,MAXNPT,3,3),
     6   	  STALL(MAXNOEL,MAXNPT,MAXSLIP),
     9      ACCGAM_ALL(MAXNOEL,MAXNPT,MAXSLIP),
     1  	  ELASJAC(MAXNOEL,MAXNPT,6,6),
	4	  RHOALL(MAXNOEL,MAXNPT,MAXSLIP),
	5	  ssdall(maxnoel,maxnpt,maxslip),
	8	  gndall(maxnoel,maxnpt,maxslip),
	7	  B0(MAXSLIP,3),C0(MAXSLIP,3),
	8	  B00(MAXSLIP,3),C00(MAXSLIP,3)  !twq
     
	real(8) QMAT(MAXNOEL,MAXNPT,3,3),
     1      FTALL(MAXNOEL,MAXNPT,3,3)

      real(8) C11,C12,C44
      real(8) GDOT0,AM,S0,H0,SS,AHARD,QL,QS
      real(8) rh0,dk,bs,ass,asu,shm,decmp_mech,yc
	real(8) IHARD,itype
	real(8) FREQTEX
      real(8) taustar   !twq
C
      INTEGER ITERK,ITERL,ITERERR
C
	real(8) RSHEAR(MAXNOEL,MAXNPT,MAXSLIP),
	4	bftt(maxnoel,maxslip)
      integer gtyp(maxnoel,maxslip),nte(maxnoel)
	
	integer NK(MAXNOEL),
	1	  iegr(MAXGRAIN,2),nb(maxgrain,maxrow),nbind(maxgrain)
	real(8) clength(maxnoel,maxslip),
	1	htrans(maxnoel,maxslip)
 	integer KINC1,NOEL1,NPT1

       integer imeg,icheck,itex,icoord,iastrain,ilstrain,iiter,idata  !twq

      End Module modname    

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1  RPL,DDSDDT,DRPLDE,DRPLDT,
     2  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED, CMNAME,
     3  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      use modname
c    
	IMPLICIT REAL*8(A-H,O-Z)

	COMMON/HARDBCC/HF0,HK1,HK2,HKP1,HKP2,HKS0
      
	CHARACTER*8 CMNAME
c
c	variables for ABAQUS UMAT
c
	DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
c
c	variables for the mechanical and material simulation
c
      dimension fold(3,3),fnew(3,3),fpinew(3,3),hardnew(maxslip),
     1 pk2new(3,3),signew(3,3),gammanew(maxslip),signewvec(6),
     2 epjac(6,6),sigtemp(3,3)

c	deformation gradient and RU decompositions	
	dimension fs(3,3),rs(3,3),us(3,3),q(3,3),qnew(3,3)
	dimension fst(3,3),as(3,3),asl(3,3),asc(3,3)
	dimension ft(3,3),fpi(3,3),fe(3,3),r(3,3),u(3,3)

      real(8),allocatable,save:: temps(:,:,:),tempn(:,:,:)

      dimension vs(maxnpt,maxslip,3),vn(maxnpt,maxslip,3)
      real(8),allocatable,save:: xl1(:,:),xl2(:,:),xl3(:,:),xc(:,:)


      real(8),allocatable,save:: ce(:,:,:)
	dimension utemps(maxgrain,maxrow,maxslip,3)
	dimension bxc(maxgrain,maxrow,3)
	dimension uvec(maxslip,3)
	dimension bc(3)
	dimension grid(3),grid2(3)
	dimension nce(maxcandi)
	dimension rrss(maxkount)
	dimension x(3),s(3),ssc0(3),ssc1(3)
	dimension hab(maxslip,maxslip),gam_inc(maxslip)
	dimension dens0(maxkount),bf(maxkount)
	dimension res_old(maxslip)

	dimension den_new(maxslip)
	dimension ssd_new(maxslip),gnd_new(maxslip)
	dimension gnd_newp(maxslip),gnd_newn(maxslip)
c-twq	integer itr,ibst,ibend,ibincr,i1,i2      
	integer itr,ibst,ibend,ibincr

      integer ngrain(maxnoel)
      real(8) gangle(MAXGRAIN,3)
      integer nelc(MAXNOEL,8)
    

!!JHK
      real(8) temptau,tempdis,vol
	real(8),allocatable,save:: bft(:)
	real(8) leng1(maxkount),
	1leng3(maxkount),leng2(maxkount),denstemp,GAML(MAXKOUNT)
	integer nnt,nn
	real(8) rrt(maxkount)
	real(8) thj,htemp,htemp1,rhot,gndt,bst
c-twq
      dimension C2CA(6,6)	
	dimension dens0g(maxkount)
	dimension templ1(maxkount),templ2(maxkount)	
     
	
	real(8),allocatable,save:: ALPHA(:,:,:),CUR_IJ(:,:,:)
	real(8),allocatable,save:: CUR_BAR(:), E0(:,:), ALLGND(:)
	real(8),allocatable,save:: TEMPBFTT(:,:)
	real(8),allocatable,save:: TEMPFTALL(:,:,:,:)
      real(8),allocatable,save:: TEMPFPITALL(:,:,:,:)
      real(8),allocatable,save:: tempssdall(:,:,:)
	real(8),allocatable,save:: tempgndall(:,:,:)  
      real(8),allocatable,save:: TEMPSTALL(:,:,:)
      real(8),allocatable,save:: TEMPTSTARTALL(:,:,:,:)         
      real(8),allocatable,save:: TEMPRHOALL(:,:,:)	
      real(8),allocatable,save:: TEMPACCGAM_ALL(:,:,:)
     	integer,allocatable,save:: INITT(:,:)         
     	real(8),allocatable,save:: VSG0(:,:,:),VNG0(:,:,:)
      integer,allocatable,save:: NSTORE(:,:)    	
      real(8),allocatable,save:: DGAM(:,:,:)
      real(8),allocatable,save:: COORD(:,:,:)
    
	
!! JHK
	if(.not.allocated(temps)) allocate(temps(maxnoel,maxslip,3))
	if(.not.allocated(tempn)) allocate(tempn(maxnoel,maxslip,3))
      if(.not.allocated(xl1)) allocate(xl1(maxnoel,maxslip))
      if(.not.allocated(xl2)) allocate(xl2(maxnoel,maxslip))
      if(.not.allocated(xl3)) allocate(xl3(maxnoel,maxslip))
      if(.not.allocated(xc)) allocate(xc(maxnoel,3))	
      if(.not.allocated(ce)) allocate(ce(maxnoel,8,3))
      if(.not.allocated(bft)) allocate(bft(maxnoel))

	if(.not.allocated(ALPHA)) allocate(ALPHA(MAXNOEL,3,3))
	if(.not.allocated(CUR_IJ)) allocate(CUR_IJ(MAXNOEL,3,3))
	if(.not.allocated(CUR_BAR)) allocate(CUR_BAR(MAXNOEL))
	if(.not.allocated(E0)) allocate(E0(MAXSLIP,3))
	if(.not.allocated(ALLGND)) allocate( ALLGND(MAXNOEL))
	if(.not.allocated(TEMPBFTT)) allocate(TEMPBFTT(MAXNOEL,MAXSLIP))
	if(.not.allocated(TEMPFTALL)) 
     1  allocate( TEMPFTALL(MAXNOEL,MAXNPT,3,3))
      if(.not.allocated(TEMPFPITALL)) 
     1  allocate( TEMPFPITALL(MAXNOEL,MAXNPT,3,3))
      if(.not.allocated(tempssdall)) 
     1 allocate( tempssdall(maxnoel,maxnpt,maxslip))
	if(.not.allocated(tempgndall)) 
     1 allocate( tempgndall(maxnoel,maxnpt,maxslip) ) 
      if(.not.allocated(TEMPSTALL)) 
     1 allocate( TEMPSTALL(MAXNOEL,MAXNPT,MAXSLIP))
      if(.not.allocated(TEMPTSTARTALL)) 
     1 allocate( TEMPTSTARTALL(MAXNOEL,MAXNPT,3,3) )        
      if(.not.allocated(TEMPRHOALL)) 
     1 allocate( TEMPRHOALL(MAXNOEL,MAXNPT,MAXSLIP)	)
      if(.not.allocated(TEMPACCGAM_ALL)) 
     1 allocate( TEMPACCGAM_ALL(MAXNOEL,MAXNPT,MAXSLIP))
     	if(.not.allocated(INITT)) allocate( INITT(10,MAXNOEL)    )     
     	if(.not.allocated(VSG0)) allocate( VSG0(MAXNOEL,MAXSLIP,3))
     	if(.not.allocated(VNG0)) allocate(VNG0(MAXNOEL,MAXSLIP,3))
      if(.not.allocated(NSTORE)) allocate( NSTORE(MAXNOEL,500)    )	
      if(.not.allocated(DGAM)) allocate( DGAM(MAXNOEL,MAXNPT,MAXSLIP))
      if(.not.allocated(COORD)) allocate( COORD(MAXNOEL,MAXNPT,3))
      
      KINC1=KINC
      NOEL1=NOEL
      NPT1=NPT
c	initialize variables
c	note that capital letter is for abaqus variables
c	lower letter is for other variables
      NSLIP=maxslip  	
      totaltime = time(2)
      ISTEP=KSTEP !STEP NUMBER FROM ABAQUS
      TIMEC=TIME(2)
c=====================================================================
      
      if(init.eq.0) then 
c
c	write to the following files      
C     
	! write messages
	open(UR7,file=Trim(AdjustL(FILEroutine))//'meg.OUT')
	! checking (temporary)
	open(UR8,file=Trim(AdjustL(FILEroutine))//'check_testa.OUT') 
	!! ip: infinity pileup, dr: drain out of free surface
	open(UR9,file=Trim(AdjustL(FILEroutine))//'iter.OUT')
	! texture data
	open(UR3,file=Trim(AdjustL(FILEroutine))//'tex_testa.OUT')
	! coordinates of each elements
	open(UR4,file=Trim(AdjustL(FILEroutine))//'coord_testa.OUT')
	! strain data
	open(UR5,file=Trim(AdjustL(FILEroutine))//'astrain.OUT')
	open(UR6,file=Trim(AdjustL(FILEroutine))//'lstrain.OUT')

	OPEN(UR10,file=Trim(AdjustL(FILEroutine))//'DATA.OUT')

c-twq  
      open(UR0,file=Trim(AdjustL(FILEroutine))//'SD.IN')
!     READS from anaqus.inp  
      call inp_info
c     OUTPUT or not?(1:Yes, 0:No)
	read(UR0,*)
	read(UR0,*)
      read(UR0,*) imeg            !imeg: output meg.txt (1:Yes, 0:No)
      read(UR0,*) icheck          !icheck: output check_testa.txt (1:Yes, 0:No)
      read(UR0,*) itex            !itex: output tex_testa.txt (1:Yes, 0:No)
      read(UR0,*) icoord          !icoord: output coord_testa.txt (1:Yes, 0:No)
      read(UR0,*) iastrain        !iastrain: output astrain.txt (1:Yes, 0:No)
      read(UR0,*) ilstrain        !ilstrain: output lstrain.txt (1:Yes, 0:No)
      read(UR0,*) iiter           !iiter: output iter.txt (1:Yes, 0:No)
      read(UR0,*) idata           !idata: output DATA.DAT (1:Yes, 0:No)   
c      
      READ(UR0,*) 
      READ(UR0,*)
      READ(UR0,*)     isimul  ! 0: only mechanical simulations 
c
      READ(UR0,*)     ihard  ! Anand:1, Nakamachi:2, density:3
      READ(UR0,*)     itype  ! 0: variable, 2: const.
      READ(UR0,*)     taustar
c     	
c	if (ihard.eq.3) then ! hardening using dislocation density
c
c     material parameters from the props of abaqus
c      
C *** READS SINGLE CRYSTAL ELASTIC STIFFNESS
      READ(UR0,*) 
      READ(UR0,*) 
      READ(UR0,*) 
      READ(UR0,*) 
      READ(UR0,*) 
      DO I=1,6
        READ(UR0,*) C2CA(I,1:6)
      ENDDO
      c11   = C2CA(1,1)  
      c12   = C2CA(1,2)   
      c44   = C2CA(4,4) 
!C *** READS PARAMETERS ABOUT CONSTITUTIVE LAW     
      READ(UR0,*)
      READ(UR0,*)     shm    ! shear modulus
      READ(UR0,*)     gdot0  ! reference shear rate
      READ(UR0,*)     am     ! rate sensitivity
      READ(UR0,*)     rh0    ! initial density
      READ(UR0,*)     bs     ! burgers vector
      READ(UR0,*)     asu,asu ! ass (Self hardening), asu (Latent hardening)
      READ(UR0,*)     dk,yc   !DK (ka:dislo.generation), yc (kb:dislo.annihilation)
c
	decmp_mech=1 ! generation decomposition

	yc=yc*bs
	rh0=rh0/MAXSLIP
	!spreading dislocation density to each slip system
       
	! initialize
c-twq	
	open(UR2,file=Trim(AdjustL(FILEroutine))//'abaqus_info.dat')   !output of abaqus.inp 
	call initial(ndi,ntens,B0,C0,VSG0,VNG0,ngrain,gangle,B00,C00)
c  
	DO IEL=1,MAXNOEL
	DO J=1,10
      INITT(J,IEL)=0     
      ENDDO
      ENDDO
	
	s0=15
			
	! element connectivity, neighboring elements
	call elcandi(NSTORE,nelc) !! element connectivity, candidates   
c-twq
      call transmission(htrans,nelc,ngrain,gangle,B00,C00)
c   
c------------------------------------------------------------------	
C     ASSIGNING GTYP VALUES FOR GRAIN BOUNDARY ELEMENTS
C
	gtyp(1:maxnoel,1:maxslip)=0 !free plane	

c	grain 1	
      read(UR2,*)      !!!GB_LIST.dat
	read(UR2,*)
	do iline=1,1
	ibkount=(iline-1)*16
	if (iline.lt.1) read (UR2,*) (nb(1,ibkount+j), j=1,16)
	if (iline.eq.1) read (UR2,*) (nb(1,ibkount+j), j=1,10)
	enddo
	nbind(1)=10
	do islip=1,nslip
	gtyp(nb(1,1:nbind(1)),islip)=2
	enddo
	
c	grain 2
	read(UR2,*)
	do iline=1,1
	ibkount=(iline-1)*16
	if (iline.lt.1) read (UR2,*) (nb(2,ibkount+j), j=1,16)
	if (iline.eq.1) read (UR2,*) (nb(2,ibkount+j), j=1,13)
	enddo
	nbind(2)=13
	do islip=1,nslip
	gtyp(nb(2,1:nbind(2)),islip)=2
	enddo

c	grain 3
	read(UR2,*)
	do iline=1,1
	ibkount=(iline-1)*16
	if (iline.lt.1) read (UR2,*) (nb(3,ibkount+j), j=1,16)
	if (iline.eq.1) read (UR2,*) (nb(3,ibkount+j), j=1,16)
	enddo
	nbind(3)=16
		do islip=1,nslip
	gtyp(nb(3,1:nbind(3)),islip)=2
	enddo

c	grain 4
	read(UR2,*)
	do iline=1,1
	ibkount=(iline-1)*16
	if (iline.lt.1) read (UR2,*) (nb(3,ibkount+j), j=1,16)
	if (iline.eq.1) read (UR2,*) (nb(3,ibkount+j), j=1,13)
	enddo
	nbind(3)=13
		do islip=1,nslip
	gtyp(nb(3,1:nbind(3)),islip)=2
	enddo
	
	statev(1) = 0.
	statev(2) = 0.

	freqtex=0.0	

      
	endif   !end if(init.eq.0)

c-------------------------------------------------------------
      strainsum=0.0
      			
      do i=1,ntens 
        strainsum = strainsum + dabs(dstran(i))
      enddo
c
      if(strainsum.lt.1.0e-15) then
        do i = 1,ntens
          do j = 1,ntens
             ddsdde(i,j) = elasjac(noel,npt,i,j)
		  enddo
		enddo 
	  return
      endif 	    	  
c
c     store the deformation gradients in ft and ftau
c
      do i=1,3
		do j=1,3
		fold(i,j)=dfgrd0(i,j)
		fnew(i,j)=dfgrd1(i,j)
		enddo
	enddo	

      IF(ISTEP.GE.2.AND.INITT(ISTEP,NOEL).LT.3)THEN

	INITT(ISTEP,NOEL)=INITT(ISTEP,NOEL)+1
	
	DO I=1,MAXNOEL
	DO J=1,MAXSLIP
	bftt(I,J)=TEMPBFTT(I,J)
	ENDDO
	ENDDO
	
      do iel=1,maxnoel
	do ipt=1,maxnpt

		do i = 1,3
		do j = 1,3
	        ftall(iel,ipt,i,j) =
     &            	tempftall(iel,ipt,i,j)
	        fpitall(iel,ipt,i,j) =
     &            	tempfpitall(iel,ipt,i,j)
	        tstartall(iel,ipt,i,j) =
     &	            temptstartall(iel,ipt,i,j)
        enddo
		enddo
	    
		do islip = 1,nslip
	       stall(iel,ipt,islip) = 
     &              tempstall(iel,ipt,islip)
             accgam_all(iel,ipt,islip) =
     &              tempaccgam_all(iel,ipt,islip)

			 rhoall(iel,ipt,islip)=
	1				temprhoall(iel,ipt,islip)
			 ssdall(iel,ipt,islip)=
	1				tempssdall(iel,ipt,islip)
			 gndall(iel,ipt,islip)=
	1				tempgndall(iel,ipt,islip)      
          enddo

	enddo
	enddo
	
	ENDIF


      if (kinc .gt. incr) then ! for new time step and kinc>1 

	if (icheck.eq.1) write(UR8,*) 'incr=', kinc
			
	
c----------------------------------------------------------------
c     STORE VARIABLES..
c----------------------------------------------------------------

	bftt(1:maxnoel,1:maxslip)=0
   
      do iel=1,maxnoel
	do ipt=1,maxnpt

		do i = 1,3
		do j = 1,3
	        ftall(iel,ipt,i,j) =
     &            	tempftall(iel,ipt,i,j)
	        fpitall(iel,ipt,i,j) =
     &            	tempfpitall(iel,ipt,i,j)
	        tstartall(iel,ipt,i,j) =
     &	            temptstartall(iel,ipt,i,j)
          enddo
		enddo
	    
		do islip = 1,nslip
	       stall(iel,ipt,islip) = 
     &              tempstall(iel,ipt,islip)
             accgam_all(iel,ipt,islip) =
     &              tempaccgam_all(iel,ipt,islip)

			 rhoall(iel,ipt,islip)=
	1				temprhoall(iel,ipt,islip)
			 ssdall(iel,ipt,islip)=
	1				tempssdall(iel,ipt,islip)
			 gndall(iel,ipt,islip)=
	1				tempgndall(iel,ipt,islip)


          enddo

	enddo
	enddo         

c----------------------------------------------------------------
      
	DO ie=1,maxnoel
		
		temps(ie,1:nslip,1:3)=0.
		tempn(ie,1:nslip,1:3)=0.
		
		do ip=1,maxnpt
			
			do i=1,3
				do j=1,3
					ft(i,j)=ftall(ie,ip,i,j)
					fpi(i,j)=fpitall(ie,ip,i,j)
				enddo
			enddo

			call mprod(ft,fpi,fe)
			call skinem(fe,r,u) ! ru decomp

			do islip=1,nslip
				do ii=1,3
					vs(ip,islip,ii)=0.
					vn(ip,islip,ii)=0.
					do jj=1,3
						vs(ip,islip,ii)=vs(ip,islip,ii)+
	1						r(ii,jj)*vsg0(ie,islip,jj)
						vn(ip,islip,ii)=vn(ip,islip,ii)+
	1						r(ii,jj)*vng0(ie,islip,jj)
					enddo
				enddo
			enddo !ISLIP

			temps(ie,1:nslip,1:3)=temps(ie,1:nslip,1:3)+
	1			vs(ip,1:nslip,1:3)
			tempn(ie,1:nslip,1:3)=tempn(ie,1:nslip,1:3)+
	1			vn(ip,1:nslip,1:3)
		enddo

		temps(ie,1:nslip,1:3)=temps(ie,1:nslip,1:3)/maxnpt !averaged
		tempn(ie,1:nslip,1:3)=tempn(ie,1:nslip,1:3)/maxnpt

	enddo

c----------------------------------------------------------------

	call getcoord(temps,tempn,xl1,xl2,xl3,xc,ce,COORD) 	
	
	if(isimul.eq.0) goto 92 !MECHANICAL SIMULATION ONLY
	

		
	do ig=1,maxgrain
		
		!retrieve grain boundary element number
	
		do irow=1,nbind(ig)
			nelem=nb(ig,irow)
			do islip=1,nslip
				stemp=0.
				do i=1,3
					stemp=stemp+
	1					temps(nelem,islip,i)*temps(nelem,islip,i)  !size
				enddo

	  	utemps(ig,irow,islip,1:3)=temps(nelem,islip,1:3)/dsqrt(stemp) ! unit
			
			enddo
			
			do i=1,3
				bxc(ig,irow,i)=xc(nelem,i)
			enddo
		enddo
		
	enddo
      
c--------------------------------------------------------------------------                                    

	do ig=1,maxgrain !For each grain 
c                       ++++++++++++++		
		if (icheck.eq.1) write(UR8,*) 'grain=', ig
	
		nbn=nbind(ig) !total number of grain boundary elements for each grain

		do ib=1,nbn !for each grain boundary element
c                      +++++++++++++++++++++++++++++++
			nelem=nb(ig,ib)

				if (icheck.eq.1) write(UR8,*) 'ib,nelem',ib,nelem

			do islip=1,nslip ! retrieve the center and unit vector				
				do i=1,3
					uvec(islip,i)=utemps(ig,ib,islip,i) ! direction
				enddo			
			enddo

			do i=1,3
				bc(i)=bxc(ig,ib,i) ! center position
			enddo	

			do islip=1,nslip	!for each slip system
c                                 +++++++++++++++++++++	
	            dgam_bc=0.
			   do ipt=1,maxnpt
				  dgam_bc=dgam_bc+dgam(nelem,ipt,islip)
			   enddo
			   
			   dgam_bc=dgam_bc/maxnpt
		   
c			   if(dabs(dgam_bc)/dtime.lt.1.e-10) then
c					goto 407
c			   endif 

c------------FIND ELEMENTS ALONG SLIP PLANE----------------------------------

			   itag=0
			   kount=1
			   xq=0.
			   dxq=xl3(nelem,islip)*1.0
			   mount=0
			   itry=1
			   ni=0
			   nnte=0
				
				do while(1)

408					continue										
					
					do i=1,3
					grid(i)=bc(i)+xq*uvec(islip,i)
					enddo

					if(kount.eq.1) then
					nte(kount)=nelem !target element

					else
					ni=nte(kount-1)  !initial guess
					endif

					if(kount.ne.1) then

						if(kount.eq.2) then
							
							ntgt=nk(ni) ! # of target elements
							nce(1:ntgt)=nstore(ni,1:ntgt) ! target element surrounding n

							i0=iegr(ig,1) !start
							i9=iegr(ig,2) !end

420							continue
							do i=1,ntgt
												
								if(nce(i).lt.i0.or.nce(i).gt.i9) then 
									ntgt=ntgt-1
									do j=i,ntgt
										nce(j)=nce(j+1)
									enddo
									goto 420
								endif
								
							enddo															
							
							call inout(ce,nce,ntgt,grid,nnte,itag)


							if(itag.eq.1.and.itry.eq.1) then ! no target, just skip
c								write(UR8,*) '1st,reversed!, itry=2'
								itry=2
								uvec(islip,1:3)=-uvec(islip,1:3)
								goto 408 
							
							endif

							if(itag.eq.1.and.itry.eq.2) then
c								write(UR8,*) 'first search fails,skip!'
								ngrid=1
c								goto 407 !option1
								goto 406 ! new
							endif
					
							if((itag.eq.0).and.(nnte.eq.nelem)) then
c							write(UR8,*) 'find but bc element,size*1.1'
								dxq=dxq*1.1
								xq=xq+dxq

c								write(UR8,*) 'ext. size by factor 1.1'
								goto 408
							endif

c							write(UR8,*) 'pass'
							

						else ! kount>=3
																								
							ntgt=nk(ni)
							nce(1:ntgt)=nstore(ni,1:ntgt)	
					
							i0=iegr(ig,1) !start
							i9=iegr(ig,2) !end

421							continue
							do i=1,ntgt
																				
								if(nce(i).lt.i0.or.nce(i).gt.i9) then 
									ntgt=ntgt-1
									do j=i,ntgt
										nce(j)=nce(j+1)
									enddo										
									goto 421
								endif								
							enddo	
													
							call inout(ce,nce,ntgt,grid,nnte,itag)										
c	lim	
							if((itag.eq.0).and.(nnte.eq.nelem)) then
								if (icheck.eq.1) write(UR8,*) 'find but bc element,size*1.1'
								dxq=dxq*1.1
								xq=xq+dxq

								if (icheck.eq.1) write(UR8,*) 'ext. size by factor 1.1'
								goto 408
							endif						

						endif ! kount=2 or 3...
							
					endif
c-----------------------------------------------------------------------
					if(itag.eq.1) then  ! no further crossing element
						ngrid=kount-1
						goto 406 ! next slip system
					
					elseif(itag.eq.0) then !found element!!
						
						if(kount.eq.1) then
						nnte=nelem
                          endif   
              if (nnte.eq.(nte(kount-1))) then
                  kount=kount-1
              endif
                          
						nte(kount)=nnte

						kount=kount+1
						xq=xq+dxq					
					
					endif

                          enddo !!end do while
                          

      
c-------------------------------------------------------------------------		
406                       continue

C	Assign GB elements and pileup elements

	do igrid=1,ngrid
		nn=nte(igrid)
		IF (igrid.eq.1) then
		gtyp(nn,islip)=2
		else
		if (gtyp(nn,islip).ne.2) then
		gtyp(nn,islip)=1		
		endif
		endif
	enddo

c	Obtain volume and size of the element
		denstemp=0
		GAML(IGRID)=0
		do igrid=1,ngrid
			nn=nte(igrid)
			vol=xl1(nn,islip)*xl3(nn,islip)
			leng2(igrid)=xl1(nn,islip) !normal to slip plane
			leng3(igrid)=xl2(nn,islip) !line direction
       		leng1(igrid)=xl3(nn,islip) !slip direction
	
			DO IPT=1,MAXNPT
			GAML(IGRID)=GAML(IGRID)+(DGAM(NN,IPT,ISLIP))
			ENDDO
	
			GAML(IGRID)=GAML(IGRID)/MAXNPT
			dens0(igrid)=RHOALL(nn,1,islip)*vol
			dens0g(igrid)=gndall(nn,1,islip)*vol			
		enddo

	nn=nte(1)
	call get_rss(islip,nn,rsst) !retrieve rss

c-----------------------------------------------------------------------------------------
      
C      WRITE(UR10,*)'INCR=',KINC
      
	CALL REDIST3(NGRID,DENS0g,BF,leng1,leng2,leng3,GAML,rsst,xc,islip)   

C-----------------------------------------------------------------------------------------

	do igrid=1,ngrid
		nn=nte(igrid)
		vol=xl1(nn,islip)*xl3(nn,islip)	
		gndall(nn,1,islip)=dens0g(igrid)/vol
	enddo

c	goto 92 ok
c--------------------------------------------------------------------------------------
417			continue

			bft(1:maxnoel)=0

			do igrid=1,ngrid
				nn=nte(igrid)
				bft(nn)=bf(igrid)
			enddo
			
			do ik=1,maxnoel
				bftt(ik,islip)=bftt(ik,islip)+bft(ik)
			enddo
			
c	goto 92 ok
c------------------------------------------------------------------------------
407		continue		

			enddo ! end islip

		enddo ! ib
	enddo ! ig

c------------------------------------------------------------------------------
c------------------------------------------------------------------
c	update and continue mechanical simulations

92	continue !! 

	if(freqtex.ge.0.) then 

	if(time(2).gt.freqtex) then
	  if (itex.eq.1) then
		write(UR3,*)
		write(UR3,*) '### new texture ###'
          write(UR3,*) 'time: ', time(2), 'inc: ', kinc
          write(UR3,*)
          write(UR3,*) 'iel, theta, phi, omega (in degree)'
      endif
      
      if (icoord.eq.1)then
		write(UR4,*) '### coordinates ###'
          write(UR4,*) 'time: ', time(2), 'inc: ', kinc
          write(UR4,*)
          write(UR4,*) 'iel,x,y,z'
          do iel=1,maxnoel
			write(UR4,7111) iel,xc(iel,1),xc(iel,2),xc(iel,3)
		enddo

		write(UR4,*)
		write(UR4,*) 'lengths'

		do iel=1,maxnoel
			write(UR4,7111) iel, xl1(iel,1), xl2(iel,1),xl3(iel,1)
			  ! check coord.	
		enddo
      endif
      
      if (iastrain.eq.1)then   
		write(UR5,*) '### almansi strain ###'		
		write(UR5,*) 'time: ', time(2), 'inc: ', kinc
	endif
		
		
		


7111		format(i5,f10.5,f10.5,f10.5,f14.2)
		
		!!! temp	
	
		do iel=1,maxnoel
		
		ft(1:3,1:3)=0.
		fpi(1:3,1:3)=0.
		q(1:3,1:3)=0.	

c		! maxnpt =1 for temp

		do ipt=1,maxnpt
c		do ipt=1,1		
			do i=1,3
			do j=1,3
				ft(i,j)=ft(i,j)+ftall(iel,ipt,i,j)
				fpi(i,j)=fpi(i,j)+fpitall(iel,ipt,i,j)
				q(i,j)=q(i,j)+qmat(iel,ipt,i,j)
			enddo
			enddo
			
		enddo

		ft(1:3,1:3)=ft(1:3,1:3)/maxnpt
		fpi(1:3,1:3)=fpi(1:3,1:3)/maxnpt ! average
		q(1:3,1:3)=q(1:3,1:3)/maxnpt

		
		! get elastic part of deformation gradient
		call mprod(ft,fpi,fs)

c
c	elastic strain components
c
	! transpose
		fst(1,1)=fs(1,1)
		fst(2,2)=fs(2,2)
		fst(3,3)=fs(3,3)

		fst(1,2)=fs(2,1)
		fst(2,1)=fs(1,2)
		fst(2,3)=fs(3,2)
		fst(3,2)=fs(2,3)
		fst(1,3)=fs(3,1)
		fst(3,1)=fs(1,3)

		call mprod(fst,fs,asc) ! c
      
		as(1,1)=0.5*(asc(1,1)-1.)
		as(2,2)=0.5*(asc(2,2)-1.)
		as(3,3)=0.5*(asc(3,3)-1.)
		as(1,2)=0.5*asc(1,2)
		as(1,3)=0.5*asc(1,3)
		as(2,1)=0.5*asc(2,1)
		as(2,3)=0.5*asc(2,3)
		as(3,1)=0.5*asc(3,1)
		as(3,2)=0.5*asc(3,2) ! almansi-hemel strain tensor

	  ! log strain
c		asl(1,1:3)=0.5*log(asc(1,1:3))	
c		asl(2,1:3)=0.5*log(asc(2,1:3))
c		asl(3,1:3)=0.5*log(asc(3,1:3))	
c		
c***NOTE*******************************************************************************************
cThis sometimes gives rise to numerical error if almansi-hemel strain tensor component is negative.   
cThe log strain seems not to be necessary. So.. I just removed the calculation.
cIf the log strain is necessary, I need to make a fine tuning.
c**************************************************************************************************

		! ru decompose
		call skinem(fs,rs,us)
		
		! get new rotation matrix
		call mprod(rs,q,qnew)

		! update texture
		call uptex(qnew,theta,phi,omega,itex_err)

		! change to degress
		call deg2rad(theta,2)
		call deg2rad(phi,2)
		call deg2rad(omega,2)
	
		if(phi.lt.0.) phi=phi+360.
		if(omega.lt.0.) omega=omega+360.

c	output rotated texture
    	if (itex.eq.1)then	
		if(itex_err.eq.1) then ! error occur
      		write(UR3,*) 'error! during update texture'
			write(UR3,*) 'skip element ', iel
		endif
		
		write(UR3,7111) iel, theta, phi, omega !!
      endif    

      if (iastrain.eq.1) then
c	write strain components
		write(UR5,1111) iel,as(1,1),as(2,2),as(3,3),as(1,2),as(2,3),
	1	as(3,1)
      endif
		enddo
      
		freqtex=freqtex+0.005
		
c     write the values if total time exceed 0, 0.01, 0.02 ....			
		
	endif

	endif

!	ok

1111	format(i4,6f8.3)

	incr  = incr + 1

	endif
c     end if above correspond to the if in page 4!!
c
	  
      CALL ZEROM(signew)

      epjac(1:NTENS,1:NTENS) = 0.0             

      iterk    = 0
      iterl    = 0
      itererr  = 0
      domppc   = 0.
      
c----------------------------------------------------
	CALL INTEG(NOEL,NPT,fold,fnew,DTIME,gammanew,
     1             fpinew,hardnew,pk2new,signew,
     2             domppc,dgmax,epjac,
     3			 den_new,ssd_new,gnd_new,gnd_newp,gnd_newn,DGAM) 
     
c     	write(UR10,*)dk

c      if(den_new(1).gt.6.e5)then
c      write(UR10,*)'INC=',KINC,'ELEMENT=',NOEL,'density=',den_new(1)
c      END IF
              
	IF (itererr .eq. 1) THEN
C
        if (imeg.eq.1)then
          WRITE(UR7,*)'ITERERR = 1:'
	    WRITE(UR7,*)'KINC = ',KINC
	    WRITE(UR7,*)'NOEL = ',NOEL,' NPT = ',NPT
          WRITE(UR7,*)'......REPEATING TIME STEP'
        endif

C
C            USE THE ELASTIC JACOBIAN AND RETURN
C
          do I = 1,NTENS
          do J = 1,NTENS
               DDSDDE(I,J) = ELASJAC(NOEL,NPT,I,J)
		  enddo
		  enddo
 
          PNEWDT = 0.75
          RETURN
	ENDIF	    

      CALL PUSHSV(signew,signewvec,2,NDI,NTENS)
C
C	UPDATE THE STRESS
C
      DO I=1,NTENS
  	      STRESS(I)=signewvec(I)
	enddo
            write(6,*) 'STRESS=', stress
C
C	UPDATE VARIABLES
C
	DO i = 1,3
	DO j = 1,3
	   TEMPFTALL(NOEL,NPT,i,j)   = fnew(i,j)
	   TEMPFPITALL(NOEL,NPT,i,j)   = fpinew(i,j)
	   TEMPTSTARTALL(NOEL,NPT,i,j) = pk2new(i,j)
	enddo
	enddo
c
      DO islip = 1,nslip
         
	   TEMPSTALL(NOEL,NPT,islip)=hardnew(islip)
         TEMPACCGAM_ALL(NOEL,NPT,islip)=
     &                                gammanew(islip)

	   TEMPRHOALL(NOEL,NPT,ISLIP)=den_new(islip)
	   tempssdall(noel,npt,islip)=ssd_new(islip)
	   tempgndall(noel,npt,islip)=gnd_new(islip)
 
	enddo
C
C	UPDATE COORDINATE SYSTEMS
C
	coord(noel,npt,1)=coords(1)
	coord(noel,npt,2)=coords(2)
	coord(noel,npt,3)=coords(3)
C
	call zerom(sigtemp)

      CALL PUSHSV(sigtemp,STRESS,1,NDI,NTENS)
      CALL EQUIVS(sigtemp,eqstress)
      IF (eqstress .EQ. 0.0) THEN
          deqpstn = 0.0
      ELSE 
          deqpstn = domppc/eqstress
      END IF
C     DOMPPC = EQUI. PLASTIC WORK (= EQUI.SIG * EQUI.EPS)
C     EQSTRESS = EQUI. STRESS (= EQUI.SIG)
C     DEQPSTRN = EQUI. PLASTIC STRAIN (=EQUI.EPS)      
C
      STATEV(1) = STATEV(1) + deqpstn
      STATEV(2) = deqpstn/DTIME   
C     STATEV(1) : EQUI. PLASTIC STRAIN
C     STATEV(2) : EQUI. PLASTIC STRAIN RATE      
C

	DO i = 1,NTENS
        DO j = 1,NTENS
          DDSDDE(I,J) = epjac(i,j)            
	  enddo
	enddo

C	STEP TIME ADJUST
C	ACCORDING TO MIT'S THESIS
C
      gammalmt= 0.03
	umeror = dgmax/gammalmt
      IF (umeror .gt. 1.) THEN
	   PNEWDT = 0.5
      ELSE IF ((umeror .gt. 0.5) .AND. (umeror .le. 1.)) THEN
         PNEWDT = 1.
      ELSE
         PNEWDT = 1.5
      END IF 

      ALLGND(NOEL)=0
c-----------------------------------------------------------------------------

      DO I=1,3
      DO J=1,3
      ALPHA(NOEL,I,J)=0
      ENDDO
      ENDDO
      
      DO I=1,3
      DO J=1,2
      CUR_IJ(NOEL,I,J)=0
      ENDDO
      ENDDO      
   
      CUR_BAR(NOEL)=0
      
      ! CALCULATION OF LINE DIRECTION VECTOR
      DO ISLIP=1,MAXSLIP
      SIZEV=dsqrt(b0(ISLIP,1)**2+b0(ISLIP,2)**2
     1 +b0(ISLIP,3)**2)
      SIZEn=dsqrt(c0(ISLIP,1)**2+c0(ISLIP,2)**2
     1 +c0(ISLIP,3)**2)
      B0(ISLIP,1:3)=b0(ISLIP,1:3)/sizeV
      C0(ISLIP,1:3)=c0(ISLIP,1:3)/sizeN	
      ENDDO
      
            
      DO ISLIP=1,MAXSLIP  
      e0(islip,1:3)=0   
      E0(ISLIP,1)=B0(islip,2)*C0(islip,3)
     1       -B0(islip,3)*C0(islip,2)
	E0(ISLIP,2)=B0(islip,3)*C0(islip,1)
	1       -B0(islip,1)*C0(islip,3)
	E0(ISLIP,3)=B0(islip,1)*C0(islip,2)
	1       -B0(islip,2)*C0(islip,1)
	size=dsqrt(E0(ISLIP,1)**2+E0(ISLIP,2)**2+E0(ISLIP,3)**2)
	E0(ISLIP,1:3)=E0(ISLIP,1:3)/size	
	ENDDO
      
      ! CALCULATION OF NYE'S TENSOR
      DO ISLIP=1,MAXSLIP
      DO IPT=1,MAXNPT
      ALPHA(NOEL,1,1)=ALPHA(NOEL,1,1)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,1)*E0(ISLIP,1)
      ALPHA(NOEL,2,2)=ALPHA(NOEL,2,2)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,2)*E0(ISLIP,2)
      ALPHA(NOEL,1,2)=ALPHA(NOEL,1,2)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,1)*E0(ISLIP,2)  
      ALPHA(NOEL,2,1)=ALPHA(NOEL,2,1)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,2)*E0(ISLIP,1)  
      ALPHA(NOEL,3,1)=ALPHA(NOEL,3,1)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,3)*E0(ISLIP,1) 
      ALPHA(NOEL,3,2)=ALPHA(NOEL,3,2)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,3)*E0(ISLIP,2)
      ALPHA(NOEL,3,3)=ALPHA(NOEL,3,3)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,3)*E0(ISLIP,3)
      ALPHA(NOEL,1,3)=ALPHA(NOEL,1,3)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,1)*E0(ISLIP,3)
      ALPHA(NOEL,2,3)=ALPHA(NOEL,2,3)+(GNDALL(NOEL,1,ISLIP))
     1       *bs*B0(ISLIP,2)*E0(ISLIP,3)
      ENDDO
      ENDDO
      
      ! AVERAGING CONSIDERING INTEGRATION POINT
      ALPHA(NOEL,1,1)=ALPHA(NOEL,1,1)
      ALPHA(NOEL,2,2)=ALPHA(NOEL,2,2)
      ALPHA(NOEL,1,2)=ALPHA(NOEL,1,2)
      ALPHA(NOEL,2,1)=ALPHA(NOEL,2,1)
      ALPHA(NOEL,3,1)=ALPHA(NOEL,3,1)
      ALPHA(NOEL,3,2)=ALPHA(NOEL,3,2)
      ALPHA(NOEL,3,3)=ALPHA(NOEL,3,3)
      ALPHA(NOEL,1,3)=ALPHA(NOEL,1,3)
      ALPHA(NOEL,2,3)=ALPHA(NOEL,2,3)
       
      ! CALCULATION OF LATTICE CURVATURE
      DO I=1,3
      DO J=1,2
      IF(I.EQ.J)THEN
      CUR_IJ(NOEL,I,J)=-ALPHA(NOEL,I,J)
     1 +(ALPHA(NOEL,1,1)+ALPHA(NOEL,2,2)+ALPHA(NOEL,3,3))/2
      ELSE
      CUR_IJ(NOEL,I,J)=-ALPHA(NOEL,I,J)
      ENDIF
      ENDDO
      ENDDO
      
      DO I=1,3
      DO J=1,2
      CUR_BAR(NOEL)=CUR_BAR(NOEL)+abs(CUR_IJ(NOEL,I,J))
      ENDDO
      ENDDO  

      DO ISLIP=1,MAXSLIP
      ALLGND(NOEL)=ALLGND(NOEL)+abs(GNDALL(NOEL,1,ISLIP))
      ENDDO  
      ALLGND(NOEL)=ALLGND(NOEL)/MAXSLIP 
      
      STATEV(11)=abs(CUR_BAR(NOEL))/6

      STATEV(4)=RHOALL(NOEL,NPT,1)    !DISLOCATIO DENSITY IN SLIP SYSTEM1

      STATEV(10)=ALLGND(NOEL)
      STATEV(12)=ALPHA(NOEL,1,1)
      STATEV(13)=ALPHA(NOEL,2,1)
      STATEV(14)=ALPHA(NOEL,3,1)
      STATEV(15)=ALPHA(NOEL,1,2)
      STATEV(16)=ALPHA(NOEL,2,2)
      STATEV(17)=ALPHA(NOEL,3,2)
      STATEV(18)=ALPHA(NOEL,1,3)
      STATEV(19)=ALPHA(NOEL,2,3)
      STATEV(20)=ALPHA(NOEL,3,3)

	
	RETURN
      END

      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++SUBROUTINES+++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-twq
      SUBROUTINE INITIAL(NDI,NTENS,VS0,VN0,VSG0,VNG0,ngrain1,gangle1,
     1    vs00,vn00) 
          use modname 
	IMPLICIT REAL*8(A-H,O-Z)	
  
!	Include 'D:\temp-wagoner102\COMMON'

      DIMENSION  rot(3,3),Qtemp(3,3),
	1 elas(3,3,3,3),elasr(6,6),ejac(6,6),rot_temp(3,3)
c-twq      real(8),allocatable:: gangle1(:,:)
      real(8) vs0(MAXSLIP,3),vn0(MAXSLIP,3)
      real(8) vs00(MAXSLIP,3),vn00(MAXSLIP,3)
!JHK	DIMENSION NG(MAXGRAIN,MAXNOEL),NUMEL(MAXGRAIN),ngrain(maxnoel)
      DIMENSION NUMEL(MAXGRAIN)
      integer,allocatable:: NG(:,:)
c      integer,allocatable:: egrain(:)      
	dimension s0_temp(maxslip),temp_s(maxslip)
	dimension ph1_temp(maxgrain),ph_temp(maxgrain),ph2_temp(maxgrain)
	DIMENSION VL0(MAXSLIP,3)
	dimension VSG0(MAXNOEL,MAXSLIP,3),VNG0(MAXNOEL,MAXSLIP,3)
      
      integer ngrain1(maxnoel)
      real(8) gangle1(MAXGRAIN,3),ph1_t,ph_t,ph2_t
	integer ig_temp
      
      allocate(NG(MAXGRAIN,MAXNOEL))
c      allocate(ngrain1(maxnoel))
      
      init=1
c
	incr=1

C-twq      PI=4.0*atan(1.0)

	CALL ONEM(ONE)
	CALL ZEROM(ZERO)

C	SET ROTATION MATRIX AS UNIT MATRIX
C	AS AN INITIAL SETTING
C
      DO iel = 1,MAXNOEL
	DO ipt  = 1,MAXNPT

        do i = 1,3
        do j = 1,3
	        rot_temp(i,j) = one(i,j)
              QMAT(iel,ipt,i,j)=ONE(i,j)
	  enddo
	  enddo

	enddo
	enddo

C
C	SET THE INITIAL ELASTIC JACOBIAN ZERO
C	
      DO iel = 1,MAXNOEL
	DO ipt  = 1,MAXNPT
		do i = 1,6
		do j = 1,6
	         ejac(i,j) = 0.0 
               ELASJAC(iel,ipt,i,j) = 0.0 
		enddo
		enddo
	enddo
	enddo

C
C	READ ORIENTATION OF CRYSTAL BY EULER ANGLES
C	AND SAVE THIS INTO THE SPECIFIC ROTATION MATRIX
C	AND CALCULATE ELASTIC MODULUS IN A GLOBAL GOORDINATE
C	SYSTEM
C
		     
      
888	continue
c-----------------------------------------------------------------
	do ilayer=1,1
		read(UR2,*)
		read(UR2,*)
		do ig=1,MAXGRAIN
			read(UR2,*) ig_temp,is_temp,ie_temp
c-twq
              iegr(ig_temp,1)=is_temp !start element
			iegr(ig_temp,2)=ie_temp !end element
c
			ngrain1(is_temp:ie_temp)=ig_temp

	    enddo

	    read(UR0,*)
		read(UR0,*)
          read(UR0,*)
		read(UR0,*)
		do ig=1,MAXGRAIN
			read(UR0,*) ig_temp,ph1_t,ph_t,ph2_t
c			ph1_temp(ig_temp)=ph1_t
c			ph_temp(ig_temp)=ph_t
c			ph2_temp(ig_temp)=ph2_t

			gangle1(ig_temp,1)=ph1_t
			gangle1(ig_temp,2)=ph_t
			gangle1(ig_temp,3)=ph2_t
	
		enddo
	enddo      
 
C
      do iel = 1,MAXNOEL	

C	retrieve the grain information
	
	id=ngrain1(iel) ! make element set for each grain
	
      do ipt = 1,MAXNPT	
 
c	 if(id.eq.1) then !! Make id as many as number of different grains

c		 ph1=ph1_temp(id)
c		 ph=ph_temp(id)
c		 ph2=ph2_temp(id)

		 ph1=gangle1(id,1)
		 ph=gangle1(id,2)
		 ph2=gangle1(id,3)

		 phi=180.-ph2
		 theta=ph
		 omega=180.-ph1

		if(phi.ge.360.) phi=phi-360.
		if(phi.lt.0.) phi=phi+360.
		if(omega.ge.360.) omega=omega-360.
		if(omega.lt.0.) omega=omega+360.
c
		  CALL DEG2RAD(theta,1)
		  CALL DEG2RAD(phi,1)
	    CALL DEG2RAD(omega,1)

	    CALL ROTMAT(theta,phi,omega,rot)

c		  write(*,*) theta,phi,omega

          do i = 1,3   
          do j = 1,3 		             
			 QMAT(iel,ipt,i,j)  = rot(i,j)
		  enddo
		  enddo
C
            CALL  ELASMAT(rot,elas)
	      CALL  EMATR(elas,elasr)
	      do i = 1,6
	      do j = 1,6
	          
			  ELASJAC(iel,ipt,i,j) = ELASJAC(iel,ipt,i,j) +
     +					    elasr(i,j)
		  enddo
		  enddo	   

	enddo
	enddo 
	
C      
C	READ SLIP SYSTEMS
C
C	GENERAL SLIP SYSTEMS
C	DIFFERENT VALUES NSLIP BETWEEN
C	BCC AND FCC
C
      read(UR0,*)
      read(UR0,*)
	read(UR0,*) nslip ! first line

c
 	do islip = 1,MAXSLIP
	      read(UR0,*)(vn00(islip,j),j=1,3),(vs00(islip,j),j=1,3)
c		  write(UR8,*)(vn0(islip,j),j=1,3),(vs0(islip,j),j=1,3)
	enddo
	
	! slip normal/slip direction
	      


	do islip = 1,NSLIP
	   
	   xxs = dSQRT(vs00(islip,1)**2+vs00(islip,2)**2+vs00(islip,3)**2)
	   xxn = dSQRT(vn00(islip,1)**2+vn00(islip,2)**2+vn00(islip,3)**2)
	   
	   do j = 1,3
	       vs0(islip,j) = vs00(islip,j)/xxs
             vn0(islip,j) = vn00(islip,j)/xxn
	   enddo          
C        Normmalization of vs0 and vn0
C
         do j=1,3
         do k=1,3
	         SMATC(islip,j,k) = vs0(islip,j)*vn0(islip,k)
		 enddo
	   enddo      

	enddo


c	calculate the slip system in global coordinate system
c	for material simulations
c	Note that this is element variables because of multicrystals
c	
	do iel=1,MAXNOEL
	do ipt=1,1 ! the same Euler angle for integration points
	
		do i=1,3
		do j=1,3
			Qtemp(i,j)=QMAT(iel,ipt,i,j)		
		enddo
		enddo

		do islip=1,nslip
		
		do i=1,3
			vsg0(iel,islip,i)=0.
			vng0(iel,islip,i)=0.
		do j=1,3

		vsg0(iel,islip,i)=vsg0(iel,islip,i)+Qtemp(i,j)*vs0(islip,j)
		vng0(iel,islip,i)=vng0(iel,islip,i)+Qtemp(i,j)*vn0(islip,j)
		
		enddo

		enddo
 
		enddo ! islip

	enddo
	enddo ! iel

C
C	SAVE INITIAL VARIABLES

      do iel = 1,MAXNOEL
      do ipt = 1,MAXNPT

            do i = 1,3
            do j = 1,3         
                 FTALL(iel,ipt,i,j)= ONE(i,j)
                 FPITALL(iel,ipt,i,j)= ONE(i,j)
		       TSTARTALL(iel,ipt,i,j)= ZERO(i,j)
			enddo
			enddo
	enddo
	enddo

      do iel = 1,MAXNOEL
      do ipt = 1,MAXNPT
         do islip = 1,NSLIP
            if(ihard.ne.3) STALL(iel,ipt,islip) = S0
c				  rhoall(iel,ipt,islip)= rh0 !! new
	   enddo
	enddo
	enddo
c
cnew	!S0
	if(ihard.eq.3) then
		
		
		if(nslip.eq.1) then
			temp_s(nslip)=rh0
		
		else
		
		do islip=1,nslip
			temp_s(islip)=0.
			do jslip=1,nslip
				
				temp1=asu
				if(jslip.eq.islip) temp1=ass
				
				temp_s(islip)=temp_s(islip)+temp1*rh0

			enddo
		enddo
		
		endif
		
c	tempres2=3.6e8
	tempres2=0
		
		do islip=1,nslip
c	  s0_temp(islip)=103+0.4*shm*bs*(temp_s(islip)+tempres2)**(1./2.)
        s0_temp(islip)=0.4*shm*bs*(temp_s(islip)+tempres2)**(1./2.)
c        s0_temp(islip)=10
		enddo
		
		do iel = 1,MAXNOEL
			do ipt = 1,MAXNPT
				do islip = 1,NSLIP
					STALL(iel,ipt,islip) = S0_temp(islip)
					rhoall(iel,ipt,islip)=rh0 !!
					ssdall(iel,ipt,islip)=rh0 
c					gndall(iel,ipt,islip)=0.5*rh0 

					gndall(iel,ipt,islip)=0
c					ssdall(iel,ipt,islip)=rh0 
c					gndall(iel,ipt,islip)=0 
				enddo
			enddo
		enddo

		s0=s0_temp(1)

	endif

	if (icheck.eq.1) write(UR8,*) 's0=',s0
C

      RETURN
      END

C--------------------------------------------------------------------------------

	subroutine getcoord(vs,vn,xl1,xl2,xl3,xc,ce,COORD1)
      use modname
	implicit real*8(a-h,o-z)
C
!	Include 'D:\temp-wagoner102\COMMON'
C
!JHK
!	dimension sh(8),gs(8,3),a(8,8),b(8),indx(8),ta(8,8)
!	dimension ce(MAXNOEL,8,3),xc(MAXNOEL,3)
!	dimension x0(3),x1(3),x2(3),x3(3),x4(3),elength(MAXNOEL,MAXSLIP)
!	dimension u(3),xl1(MAXNOEL,MAXSLIP),xl2(MAXNOEL,MAXSLIP),
!	1		  xl3(MAXNOEL,MAXSLIP)
!	dimension vs(MAXNOEL,MAXSLIP,3),vn(MAXNOEL,MAXSLIP,3)
!	dimension va(3),vb(3)
!      DIMENSION COORD1(MAXNOEL,MAXNPT,3)	

      dimension sh(8),gs(8,3),a(8,8),b(8),indx(8),ta(8,8)
	dimension ce(MAXNOEL,8,3),xc(MAXNOEL,3)
	dimension x0(3),x1(3),x2(3),x3(3),x4(3)
	real(8),allocatable:: elength(:,:)
	dimension u(3),xl1(MAXNOEL,MAXSLIP),xl2(MAXNOEL,MAXSLIP),
	1		  xl3(MAXNOEL,MAXSLIP)
	dimension vs(MAXNOEL,MAXSLIP,3),vn(MAXNOEL,MAXSLIP,3)
	dimension va(3),vb(3)
      DIMENSION COORD1(MAXNOEL,MAXNPT,3)	

	allocate(elength(MAXNOEL,MAXSLIP))
C
C	gauss points
C	
	gs(1,1)=-1./dsqrt(3.d0) ; gs(1,2)=-1./dsqrt(3.d0) ; 
	gs(1,3)=-1./dsqrt(3.d0)
	gs(2,1)= 1./dsqrt(3.d0) ; gs(2,2)=-1./dsqrt(3.d0) ; 
	gs(2,3)=-1./dsqrt(3.d0)
	gs(3,1)=-1./dsqrt(3.d0) ; gs(3,2)= 1./dsqrt(3.d0) ; 
	gs(3,3)=-1./dsqrt(3.d0)
	gs(4,1)= 1./dsqrt(3.d0) ; gs(4,2)= 1./dsqrt(3.d0) ; 
	gs(4,3)=-1./dsqrt(3.d0)
	gs(5,1)=-1./dsqrt(3.d0) ; gs(5,2)=-1./dsqrt(3.d0) ; 
	gs(5,3)= 1./dsqrt(3.d0)
	gs(6,1)= 1./dsqrt(3.d0) ; gs(6,2)=-1./dsqrt(3.d0) ; 
	gs(6,3)= 1./dsqrt(3.d0)
	gs(7,1)=-1./dsqrt(3.d0) ; gs(7,2)= 1./dsqrt(3.d0) ; 
	gs(7,3)= 1./dsqrt(3.d0)
	gs(8,1)= 1./dsqrt(3.d0) ; gs(8,2)= 1./dsqrt(3.d0) ; 
	gs(8,3)= 1./dsqrt(3.d0)
c
c	make array (8 x 8)
c
	a(1:8,1:8)=0.
c
	do i=1,8
	
		xi=gs(i,1)
		eta=gs(i,2)
		chi=gs(i,3)

		call shapef(xi,eta,chi,sh)

		do j=1,8

			a(i,j)=sh(j)

		enddo

	enddo

c
c	coord=a^-1*coordIntp.
c
		
	b(1:8)=0.
c
	do i=1,3 ! for x,y,z

		do iel=1,MAXNOEL
	
			ta(1:8,1:8)=a(1:8,1:8)
			
			do ipt=1,8
				
	 			b(ipt)=coord1(iel,ipt,i)
			
			enddo
	
			call ludcmp(ta,8,8,indx,d)
			call lubksb(ta,8,8,indx,b) ! b is solution
	

			do j=1,8
				ce(iel,j,i)=b(j) ! element coordinates, 1~8 :element connectivity
			enddo

		enddo

	enddo
	 
c
c	calculate center coordinate
c	
	
	do ie=1,MAXNOEL

		xi=0.
		eta=0.
		chi=0.

		call shapef(xi,eta,chi,sh)

		do j=1,3
			xc(ie,j)=0.
			
			do inode=1,8
			
				xc(ie,j)=xc(ie,j)+sh(inode)*ce(ie,inode,j)

			enddo

		enddo

	enddo
c

C
C	calculate intersection point FOR EACH SLIP SYSTEM
C	NOTE: THE VN,VS ARE NORMAL VECTOR IN THE GLOBAL COORDINATE
C		  SYSTEM	
C
	icase=1

101	continue
c
	do ie=1,MAXNOEL ! for each element
	
	do  islip=1,nslip

	if(icase.eq.1) then
		
		u(1:3)=vn(ie,islip,1:3)

	else if(icase.eq.2) then

		va(1:3)=vn(ie,islip,1:3)
		vb(1:3)=vs(ie,islip,1:3)
		
		call crospd(va,vb,u)

	else if(icase.eq.3) then

		u(1:3)=vs(ie,islip,1:3)

	endif

	xuu=dsqrt(u(1)**2+u(2)**2+u(3)**2)

	u(1:3)=u(1:3)/xuu	
	
	elength(ie,islip)=0.
	icount=0

	do i=1,6 ! six planes in single cubic

		x0(1:3)=xc(ie,1:3) ! center

		if(i.eq.1) then
		
			x1(1:3)=ce(ie,1,1:3)
			x2(1:3)=ce(ie,5,1:3)
			x3(1:3)=ce(ie,6,1:3)
			x4(1:3)=ce(ie,2,1:3)

		elseif(i.eq.2) then

			x1(1:3)=ce(ie,4,1:3)
			x2(1:3)=ce(ie,3,1:3)
			x3(1:3)=ce(ie,7,1:3)
			x4(1:3)=ce(ie,8,1:3)


		elseif(i.eq.3) then

			x1(1:3)=ce(ie,6,1:3)
			x2(1:3)=ce(ie,2,1:3)
			x3(1:3)=ce(ie,3,1:3)
			x4(1:3)=ce(ie,7,1:3)

		elseif(i.eq.4) then

			x1(1:3)=ce(ie,5,1:3)
			x2(1:3)=ce(ie,1,1:3)
			x3(1:3)=ce(ie,4,1:3)
			x4(1:3)=ce(ie,8,1:3)

		elseif(i.eq.5) then

			x1(1:3)=ce(ie,1,1:3)
			x2(1:3)=ce(ie,2,1:3)
			x3(1:3)=ce(ie,3,1:3)
			x4(1:3)=ce(ie,4,1:3)

		elseif(i.eq.6) then

			x1(1:3)=ce(ie,5,1:3)
			x2(1:3)=ce(ie,6,1:3)
			x3(1:3)=ce(ie,7,1:3)
			x4(1:3)=ce(ie,8,1:3)


		endif
			
		itag=0
				
		call getintersect(i,u,x0,x1,x2,x3,x4,halfL,itag)

c		write(UR10,*) 'itag',itag

		if(itag.ge.1) then
		
			icount=icount+1
			elength(ie,islip)=elength(ie,islip)+dabs(halfL)

c			write(UR10,*) 'icount, halfL', icount,halfL
		
		endif

		if(icount.eq.2) goto 100
				
		if(i.eq.6.and.icount.lt.2) then

			if (icheck.eq.1) write(UR8,*) 'Error in obtaining Length'
			stop

		endif

	enddo ! i

100	continue ! get length for this element


	enddo ! islip

	enddo ! element


	if(icase.eq.1) then

		xl1(1:MAXNOEL,1:NSLIP)=elength(1:MAXNOEL,1:NSLIP)

		icase=icase+1

		goto 101

	elseif(icase.eq.2) then

		xl2(1:MAXNOEL,1:NSLIP)=elength(1:MAXNOEL,1:NSLIP)
		
		icase=icase+1

		goto 101

	else

		xl3(1:MAXNOEL,1:NSLIP)=elength(1:MAXNOEL,1:NSLIP)

	endif
	
	do iel=1,maxnoel
		do islip=1,nslip
			clength(iel,islip)=xl3(iel,islip)
		enddo
	enddo

c	do ie=1,maxnoel
c		do islip=1,nslip
c			xl1(ie,islip)=0.0001
c			xl2(ie,islip)=0.0001
c		enddo
c	enddo
C
C
C
	return
	end
c
c--------------------------------------------------------------------------
c
	subroutine getintersect(iface,u,x0,x1,x2,x3,x4,halfL,itag)

	implicit real*8(a-h,o-z)
c
	dimension x0(3),x1(3),x2(3),x3(3),x4(3)
	dimension xnorm(3),t(3),u(3),tempv1(3)
	dimension x(3),va(3),vb(3),xt(3)
c

c	calculate normal for the lst triangle

	va(1:3)=x2(1:3)-x1(1:3)
	vb(1:3)=x4(1:3)-x1(1:3)
	
	call crospd(va,vb,xnorm)

c	find intersection

	tempv1(1:3)=x1(1:3)-x0(1:3)
	
	call dotpd(tempv1,xnorm,temps1)			
	
	call dotpd(u,xnorm,temps2)

	if(abs(temps2).le.1.e-10) return
		
	xu=temps1/temps2

	x(1:3)=x0(1:3)+xu*u(1:3) ! intersection

c	write(UR10,*) 'inter', x

c
c	check if it is inner or not
c
	
	xt(1:3)=x(1:3)-x1(1:3)

	! 1-2
	
	t1=va(1)
	t2=va(2)
	s1=vb(1)
	s2=vb(2)

	xx1=xt(1)
	xx2=xt(2)
	
	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) goto 200
	
	! 2-3	

	t1=va(2)
	t2=va(3)
	s1=vb(2)
	s2=vb(3)

	xx1=xt(2)
	xx2=xt(3)

	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) goto 200
	
	! 3-1
		
	t1=va(3)
	t2=va(1)
	s1=vb(3)
	s2=vb(1)

	xx1=xt(3)
	xx2=xt(1)

	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) then
	
		goto 200

	else

		if (idata.eq.1) write(UR10,*) 'Error'
		stop
	
	endif

200	continue

	p=((s2-t2)*xx1+(t1-s1)*xx2)/xjacob
	q=(-t2*xx1+t1*xx2)/xjacob


c	write(UR10,*) 'p,q',p,q

	if((p.ge.-1.e-6.and. p.le.(1.+1.e-6)).and.(q.le.(p+1.e-6))
	1		.and. (q.ge.-1.e-6).and.(q.le.(1.+1.e-6))) then

		itag=itag+1

		halfL=xu

		return

	endif

	continue ! neighbor triangle

c	calculate normal for the lst triangle

	va(1:3)=x3(1:3)-x2(1:3)
	vb(1:3)=x4(1:3)-x2(1:3)
	
	call crospd(va,vb,xnorm)
	
c	find intersection

	tempv1(1:3)=x2(1:3)-x0(1:3)
	
	call dotpd(tempv1,xnorm,temps1)			
	
	call dotpd(u,xnorm,temps2)

	if(abs(temps2).le.1.e-10) return

	xu=temps1/temps2

	x(1:3)=x0(1:3)+xu*u(1:3) ! intersection

c	write(UR10,*) '*inter', x


c
c	check if it is inner or not
c
	xt(1:3)=x(1:3)-x2(1:3)
	
	t1=va(1)
	t2=va(2)
	s1=vb(1)
	s2=vb(2)

	xx1=xt(1)
	xx2=xt(2)
	
	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) goto 300
	
	! 2-3	

	t1=va(2)
	t2=va(3)
	s1=vb(2)
	s2=vb(3)

	xx1=xt(2)
	xx2=xt(3)

	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) goto 300
	
	! 3-1
		
	t1=va(3)
	t2=va(1)
	s1=vb(3)
	s2=vb(1)

	xx1=xt(3)
	xx2=xt(1)

	xjacob=t1*(s2-t2)-t2*(s1-t1)
	
	if(abs(xjacob).gt.1.e-10) then
	
		goto 300

	else

		if (idata.eq.1) write(UR10,*) 'Error'
		stop
	
	endif


300	continue

	p=((s2-t2)*xx1+(t1-s1)*xx2)/xjacob
	q=(-t2*xx1+t1*xx2)/xjacob

c	write(UR10,*) '*p,q',p,q
	
	if((p.ge.-1.e-6.and. p.le.(1.+1.e-6)).and.(q.le.(p+1.e-6))
	1		.and. (q.ge.-1.e-6).and.(q.le.(1.+1.e-6))) then

		itag=itag+1

		halfL=xu

		return

	endif


	return
	end
c
C------------------------------------------------------------------------

	subroutine dotpd(va,vb,xdot)	

	implicit real*8(a-h,o-z)
	
	dimension va(3),vb(3)

	xdot=va(1)*vb(1)+va(2)*vb(2)+va(3)*vb(3)

	return
	end	
c
c
	subroutine crospd(va,vb,xnorm)

	implicit real*8(a-h,o-z)

	dimension va(3),vb(3),xnorm(3)

	xnorm(1)=va(2)*vb(3)-va(3)*vb(2)
	xnorm(2)=va(3)*vb(1)-va(1)*vb(3)
	xnorm(3)=va(1)*vb(2)-va(2)*vb(1)

	size=dsqrt(xnorm(1)**2+xnorm(2)**2+xnorm(3)**2)

	xnorm(1:3)=xnorm(1:3)/size
	
	
	return
	end
c
C------------------------------------------------------------------------
c
	subroutine shapef(xi,eta,chi,sh)
	
	implicit real*8(a-h,o-z)

	dimension sh(8)


	sh(1)=1./8.*(1.-xi)*(1.-eta)*(1.-chi)
	sh(2)=1./8.*(1.+xi)*(1.-eta)*(1.-chi)
	sh(3)=1./8.*(1.+xi)*(1.+eta)*(1.-chi)
	sh(4)=1./8.*(1.-xi)*(1.+eta)*(1.-chi)
	sh(5)=1./8.*(1.-xi)*(1.-eta)*(1.+chi)
	sh(6)=1./8.*(1.+xi)*(1.-eta)*(1.+chi)
	sh(7)=1./8.*(1.+xi)*(1.+eta)*(1.+chi)
	sh(8)=1./8.*(1.-xi)*(1.+eta)*(1.+chi)


	return
	end

C---------------------------------------------------------------------

       SUBROUTINE INTEG(NOEL,NPT,dfgold,dfgnew,DTIME,gammanew,
     &          fpinew,gnew,pk2new,signew,DOMPPC,DGMAX,WJMATGR,
     &			RHONEW,ssdnew,gndnew,gndnewp,gndnewn,DGAM1)
       use modname    	
	IMPLICIT REAL*8(A-H,O-Z)     
 

!	Include 'D:\temp-wagoner102\COMMON'
C
	DIMENSION   dfgnew(3,3),fpinew(3,3),pk2new(3,3),
	1			signew(3,3),TTAU(3,3),gnew(MAXSLIP),
	2			gammanew(MAXSLIP),FPIT(3,3),TSTART(3,3),
     3			TSTARTV(6),ST(MAXSLIP),Q(3,3),QT(3,3),
     4			SMATG(MAXSLIP,3,3),AMAT(3,3),
	5			TSTARTR(3,3),TSTARTRV(6),SALPHA(3,3),
     6			PALPHA(3,3),PALPHAV(6),PMAT(MAXSLIP,6),
	7			CALPHA(3,3),CALPHAV(6),CMAT(MAXSLIP,6),
	8			TSTARTAU(3,3),TSTARTAUV(6),STAU(MAXSLIP),
	9			RSS(MAXSLIP),DGAMMA(MAXSLIP),FPITAU(3,3),
     1			FSTAU(3,3),AUX1(3,3),AUX2(3,3),AUX3(3,3)
     	
	DIMENSION  FSTAUINV(3,3),dfgold(3,3),FTTAU(3,3),RTTAU(3,3),
	1			UTTAU(3,3),FSTART(3,3),CJMAT(3,3,3,3),
     2			ELMAT(3,3,3,3),DJMAT(3,3,3,3),DJMATR(6,6),
	3			GJMAT(MAXSLIP,3,3,3,3),AJJMAT(MAXSLIP,3,3,3,3),
	4			AJJMATR(MAXSLIP,6,6),BJMAT(MAXSLIP,3,3),
     5			BJMATR(MAXSLIP,6),AKJMATR(6,6),AKJMATRINV(6,6),
	6			AMJMATR(6,6),QJMATR(6,6),QJMAT(3,3,3,3),
	7			RJMAT(MAXSLIP,3,3),SJMAT(3,3,3,3),WJMAT(3,3,3,3),
     8			WJMATR(6,6),WJMATGR(6,6),INDX(6),
	9			AUXMAT1(3,3,3,3),AUXMAT1R(6,6),AUXMAT2(3,3,3,3),
	1			AUXV1(6),GAMMAT(MAXSLIP)

	dimension   RHOT(MAXSLIP),RHOTAU(MAXSLIP),RHONEW(MAXSLIP)
	dimension   ssdt(maxslip), ssdtau(maxslip),ssdnew(maxslip)
	dimension   gndt(maxslip), gndtau(maxslip),gndnew(maxslip)
	dimension   gndtp(maxslip), gndtaup(maxslip),gndnewp(maxslip)
	dimension   gndtn(maxslip), gndtaun(maxslip),gndnewn(maxslip)
	DIMENSION DGAM1(MAXNOEL,MAXNPT,MAXSLIP)

C	INITIALIZE
	DO I = 1,3
	DO J = 1,3
	   pk2new(I,J) = 0.
	   fpinew(I,J)   = 0.
	ENDDO
	ENDDO

	DO ISLIP = 1,NSLIP
	   gnew(ISLIP)=0.
	   gammanew(ISLIP)=0.
	ENDDO
C
	CALL ZEROM(signew)		
      ITERKMAX = 0
	ITERLMAX = 0
      DOMPPC = 0.0
C
C	ACCUMULATED GAMMA IN THE PREVIOUS STEP
C
	DO ISLIP=1,NSLIP

		GAMMAT(ISLIP)=ACCGAM_ALL(NOEL,NPT,ISLIP)

	ENDDO

	DO J = 1,3
      DO K = 1,3
         FPIT(J,K)   = FPITALL(NOEL,NPT,J,K)
	   TSTART(J,K) = TSTARTALL(NOEL,NPT,J,K)
	ENDDO
	ENDDO
		
	CALL SMATVEC(TSTART,TSTARTV)

      DO J = 1,NSLIP
         ST(J) = STALL(NOEL,NPT,J)
	ENDDO

cnew
	do j=1,nslip
	   RHOT(j)=RHOALL(NOEL,NPT,J)
	   ssdt(j)=ssdall(noel,npt,j)
	   gndt(j)=gndall(noel,npt,j)

	enddo
cnew

      DO J = 1,3
      DO K = 1,3
         Q(J,K) = QMAT(NOEL,NPT,J,K)
	ENDDO
	ENDDO

      CALL MTRANS(Q,QT)	
C
C	CALCULATE MATRIX A IN THE THESIS
C
	CALL MPROD(dfgnew,FPIT,AUX1)
	CALL MTRANS(AUX1,AUX2)
	CALL MPROD(AUX2,AUX1,AMAT)

C	CALCULATE TRIAL STRESS
	
	DO J = 1,3
	DO K = 1,3
	   AUX1(J,K) = 0.5*(AMAT(J,K) - ONE(J,K))
	ENDDO
	ENDDO
	
	CALL  ELASTOP(AUX1,Q,QT,TSTARTR)
	CALL  SMATVEC(TSTARTR,TSTARTRV)

	DO ISLIP = 1,NSLIP
C
	   DO I=1,3
	   DO J=1,3
		  SALPHA(I,J) = SMATC(ISLIP,I,J)
		  SMATG(ISLIP,I,J) = 0.0
	   ENDDO
	   ENDDO
C
C	CONVERT S INTO THE GLOBAL COORDINATE SYSTEM
C
         CALL MPROD(SALPHA,QT,AUX1)
         CALL MPROD(Q,AUX1,SALPHA)
	   DO I = 1,3
	   DO J = 1,3		   
		  SMATG(ISLIP,I,J) =SALPHA(I,J)
	   ENDDO
	   ENDDO		   
C
	   DO I = 1,3
	   DO J = 1,3
            PALPHA(I,J) = 0.5*(SALPHA(I,J)+SALPHA(J,I))
	   ENDDO
	   ENDDO
         
	   CALL SMATVEC(PALPHA,PALPHAV)

C
	   PALPHAV(4) = 2.0*PALPHAV(4)
	   PALPHAV(5) = 2.0*PALPHAV(5)
	   PALPHAV(6) = 2.0*PALPHAV(6)
	 
	   DO I=1,6
	      PMAT(ISLIP,I) = PALPHAV(I)
	   ENDDO

         CALL MPROD(AMAT,SALPHA,AUX1)
	   
	   DO I = 1,3
	   DO J = 1,3
		  AUX2(I,J) = 0.5*(AUX1(I,J) + AUX1(J,I))
	   ENDDO
	   ENDDO
	   
	   CALL ELASTOP(AUX2,Q,QT,CALPHA)
	   CALL SMATVEC(CALPHA,CALPHAV)
	   DO I=1,6
	      CMAT(ISLIP,I) = CALPHAV(I)
	   ENDDO

	ENDDO
C
	DO I = 1,6
	    TSTARTAUV(I) = TSTARTV(I)
	ENDDO

      DO I = 1,NSLIP
	   STAU(I) =ST(I)
	   DGAMMA(I) = 0.0
	ENDDO

      ITERK = 0
	ITERL = 0
	ITERERR = 0

	CALL SOLVCRYS(GAMMAT,TSTARTV,ST,TSTARTRV,PMAT,CMAT,DTIME,
	1 TSTARTAUV,STAU,DGAMMA,DGMAX,RSS,DOMPSC,ITERK,ITERL,ITERERR,
     2 RHOT,RHOTAU,ssdt,ssdtau,gndt,gndtau,NOEL)
     
c	
      IF (ITERERR .EQ. 1) THEN
          if (imeg.eq.1) WRITE(UR7,*)'ITERERR = 1! IN 
     1            INTEG AFTER CALL TO SOLVCRYS'
          RETURN
      END IF

      CALL VECSMAT(TSTARTAUV,TSTARTAU)
	DOMPPC = DOMPPC + DOMPSC
      ITERKMAX = MAX0(ITERK,ITERKMAX)
	ITERLMAX = MAX0(ITERL,ITERLMAX)

C	   
	DO I =1,3
	DO J =1,3
	   AUX1(I,J)= 0.0
		 
	DO ISLIP=1,NSLIP
	   AUX1(I,J)= AUX1(I,J) + 
     &   DGAMMA(ISLIP)*SMATG(ISLIP,I,J)                   
	ENDDO

         AUX2(I,J) = ONE(I,J) - AUX1(I,J)

	ENDDO
	ENDDO   

	CALL MPROD(FPIT,AUX2,FPITAU)
	CALL MDET(FPITAU,DET)
	CDET = DET**(1./3.)
	   
	DO I =1,3
	DO J=1,3
	   FPITAU(I,J) = FPITAU(I,J)/CDET
	ENDDO
	ENDDO
C
C	FSTAU 
C
	CALL MPROD(dfgnew,FPITAU,FSTAU)   
C
C	CAUCHY STRESS
C
	CALL MTRANS(FSTAU,AUX1)
	CALL MPROD(TSTARTAU,AUX1,AUX2)
	CALL MPROD(FSTAU,AUX2,AUX3)
	CALL MDET(FSTAU,DET)
C
	DO I = 1,3
	DO J = 1,3
	   TTAU(I,J) = AUX3(I,J)/DET
	ENDDO
	ENDDO	

C
	DO I = 1,3
	DO J = 1,3
	   SIGNEW(I,J) = SIGNEW(I,J) + TTAU(I,J)
	ENDDO
	ENDDO

C	UPDATE
C
	DO I = 1,3
	DO J = 1,3
	   PK2NEW(I,J) = TSTARTAU(I,J)
	   FPINEW(I,J) = FPITAU(I,J)
	ENDDO
	ENDDO


	DO ISLIP = 1,NSLIP
	  GNEW(ISLIP)=STAU(ISLIP)
cnew
	  RHONEW(ISLIP)=RHOTAU(ISLIP) !
	  ssdnew(islip)=ssdtau(islip)
	  gndnew(islip)=gndtau(islip)

cnew
	  GAMMANEW(ISLIP)= 
     &	     ACCGAM_ALL(NOEL,NPT,ISLIP)
     &       + DABS(DGAMMA(ISLIP))
	ENDDO
C
C	STORE DGAMMA AND RESOLVED SHEAR STRESS OF EACH SLIP SYSTEM
C	STORE AS A EXTERNAL VARIABLES(COMMON)
C	FOR MATERIAL SIMULATIONS
C
	DO ISLIP=1,NSLIP

		DGAM1(NOEL,NPT,ISLIP)=DGAMMA(ISLIP)

		RSHEAR(NOEL,NPT,ISLIP)=RSS(ISLIP)
		
C     UPDATING VALUES FOR STATE VARIABLES FOR UMAT VALIDATION
 		

        
	ENDDO
C
C
C----------------------------------------------------------------
C                 CALCULATE THE JACOBIAN
C----------------------------------------------------------------
C
C                     STEP J1. CALCULATE THE C_IJKL MATRIX
C
	      CALL M3INV(dfgold,AUX1)	      
	      CALL MPROD(dfgnew,AUX1,FTTAU)
	      CALL SKINEM(FTTAU,RTTAU,UTTAU)
	      CALL MPROD(dfgold,FPIT,FSTART)
	      
            DO 150 I = 1,3
	      DO 150 J = 1,3
	      DO 150 K = 1,3
	      DO 150 L = 1,3
	        SUM1 = 0.0
		      SUM2 = 0.0
		    DO 148 M = 1,3
		      SUM1 = SUM1 + UTTAU(L,M)  * FSTART(M,J)
		      SUM2 = SUM2 + FSTART(M,I) * UTTAU(M,K)
148         CONTINUE	      
	            CJMAT(I,J,K,L) = FSTART(K,I) * SUM1 +
     &                               SUM2 * FSTART(L,J)
150         CONTINUE	    
C   
C                  STEP J2. CALCULATE THE L_IJKL MATRIX IN GLOBAL
C                           CORDINATES
C            
            CALL  ELASMAT(Q,ELMAT)
C	      
C                  STEP J3. CALCULATE THE D_IJKL MATRIX
C
            DO 160 I = 1,3
	      DO 160 J = 1,3
	      DO 160 K = 1,3
	      DO 160 L = 1,3
	         DJMAT(I,J,K,L) = 0.0
		    DO 158 M  = 1,3
		    DO 158 N  = 1,3
		       DJMAT(I,J,K,L) = DJMAT(I,J,K,L) +
     &                     0.5*ELMAT(I,J,M,N)*CJMAT(M,N,K,L)
158         CONTINUE
160         CONTINUE
            CALL EMATR(DJMAT,DJMATR)     
C
C              STEP J4. CALCULATE THE NLIP MATRICES G_MNKL^{\ALPHA}
C
            DO 170 ISLIP = 1,NSLIP
	      DO 168 M = 1,3
		    DO 168 N = 1,3
		    DO 168 K = 1,3
		    DO 168 L = 1,3
		       GJMAT(ISLIP,M,N,K,L) = 0.0
		    DO 166 IP = 1,3
		       GJMAT(ISLIP,M,N,K,L) = GJMAT(ISLIP,M,N,K,L) +
     &                          CJMAT(M,IP,K,L)*SMATG(ISLIP,IP,N) +
     &                          SMATG(ISLIP,IP,M)*CJMAT(IP,N,K,L)
166         CONTINUE
168         CONTINUE
170         CONTINUE
C  
C              STEP J5. CALCULATE THE NSLIP MATRICES J_IJKL^{\ALPHA}
C            
	      DO 190 ISLIP = 1,NSLIP
            DO 180 I = 1,3
	      DO 180 J = 1,3
	      DO 180 K = 1,3
	      DO 180 L = 1,3
	         AJJMAT(ISLIP,I,J,K,L) = 0.0
		    DO 175 M  = 1,3
		    DO 175 N  = 1,3
		       AJJMAT(ISLIP,I,J,K,L) = AJJMAT(ISLIP,I,J,K,L) +
     &                  0.5*ELMAT(I,J,M,N)*GJMAT(ISLIP,M,N,K,L)
175         CONTINUE
               AUXMAT1(I,J,K,L) = AJJMAT(ISLIP,I,J,K,L)
180         CONTINUE
            CALL EMATR(AUXMAT1,AUXMAT1R)
		    DO 185 I = 1,6
		    DO 185 J = 1,6
		       AJJMATR(ISLIP,I,J) = AUXMAT1R(I,J)
185         CONTINUE		   
190         CONTINUE
C
C              STEP J6. CALCULATE THE NSLIP MATRICES B_IJ^{\ALPHA}
C
            DO 210 ISLIP = 1,NSLIP
	      DO 200 I     = 1,3
	   	DO 200 J     = 1,3
	         BJMAT(ISLIP,I,J) = 
     &              0.5*(SMATG(ISLIP,I,J) +SMATG(ISLIP,J,I))*
     &             GDOT0*DTIME* DABS(RSS(ISLIP))**(1.D0/AM - 1.)/
     & 		       (AM*DABS(STAU(ISLIP))**(1./AM))
		       AUX1(I,J) = BJMAT(ISLIP,I,J)
200         CONTINUE
C
            CALL SMATVEC(AUX1,AUXV1)
C
		    DO 208  I=1,6
		       BJMATR(ISLIP,I) = AUXV1(I)
208         CONTINUE
210         CONTINUE
C
C                STEP J7. CALCULATE THE  REDUCED MATRIX K_IJ
C
            DO 216 I = 1,6
	      DO 216 J = 1,6
                     AKJMATR(I,J) = 0.0
		    DO 214 ISLIP = 1,NSLIP
		       AKJMATR(I,J) = AKJMATR(I,J) +
     &                CMAT(ISLIP,I)*BJMATR(ISLIP,J)
214         CONTINUE
              IF(I .EQ. J) THEN
		       AKJMATR(I,J) = 1. + AKJMATR(I,J)
               END IF
216         CONTINUE
C
C                STEP J8. CALCULATE THE INVERSE OF THE MATRIX K_IJ
C
	      CALL MATINV(AKJMATR,6,6,INDX,AKJMATRINV)
C
C                STEP J9. CALCULATE THE REDUCED MATRIX M_IJ
C
            DO 222 I = 1,6
	      DO 222 J = 1,6
		       AMJMATR(I,J) = 0.0
		    DO 220 ISLIP= 1, NSLIP
		       AMJMATR(I,J) = AMJMATR(I,J) +
     &                     DGAMMA(ISLIP)*AJJMATR(ISLIP,I,J)
220         CONTINUE
		       AMJMATR(I,J) =  DJMATR(I,J) - AMJMATR(I,J)
222         CONTINUE
C
C                STEP J10. CALCULATE THE  MATRIX Q_IJKL
C     
            DO 226 I = 1,6
	      DO 226 J = 1,6
	         QJMATR(I,J) = 0.0
		    DO 224 M  = 1,6
               QJMATR(I,J) = QJMATR(I,J)+ AKJMATRINV(I,M)*
     &                               AMJMATR(M,J)
224         CONTINUE
226         CONTINUE
            CALL EMATREC(QJMATR,QJMAT)
C
C               STEP J11. CALCULATE THE MATRICES R_IJ^{\ALPHA}
C                 
            DO  228 ISLIP = 1, NSLIP
		    DO  228     I = 1,3
		    DO  228     J = 1,3
		        RJMAT(ISLIP,I,J) = 0.0
			DO 227 K = 1,3
			DO 227 L = 1,3
			    RJMAT(ISLIP,I,J) = RJMAT(ISLIP,I,J) +
     &                       BJMAT(ISLIP,K,L)*QJMAT(K,L,I,J)
227         CONTINUE
228                CONTINUE
C
C              STEP J12. CALCULATE THE MATRIX S_IJKL
C 
	       DO 286 IP = 1,3
	       DO 286  J = 1,3
	       DO 286  K = 1,3
	       DO 286  L = 1,3
	          AUX1(IP,J)        = 0.0
		        AUXMAT1(IP,J,K,L) = 0.0
		     DO 284 ISLIP = 1,NSLIP
		        AUX1(IP,J) = AUX1(IP,J) 
     &	                    + DGAMMA(ISLIP)*SMATG(ISLIP,IP,J)
                AUXMAT1(IP,J,K,L) = AUXMAT1(IP,J,K,L)
     &                      + SMATG(ISLIP,IP,J)*RJMAT(ISLIP,K,L)
284          CONTINUE		  
286          CONTINUE
C		    
	       DO 290 L = 1,3
	       DO 290 J = 1,3
	           AUX2(L,J) = 0.0
		     DO 288 IP = 1,3
  	           AUX2(L,J) = AUX2(L,J) +
     &                  FSTART(L,IP)*AUX1(IP,J)	      
288          CONTINUE		      
290          CONTINUE		    
C
	       DO 294 N = 1,3
	       DO 294 J = 1,3
	       DO 294 K = 1,3
	       DO 294 L = 1,3
	          AUXMAT2(N,J,K,L) = 0.0
	          DO 292  IP = 1,3		     
	            AUXMAT2(N,J,K,L) =  AUXMAT2(N,J,K,L) +
     &                      FSTART(N,IP)*AUXMAT1(IP,J,K,L)
292               CONTINUE
294            CONTINUE	
C
	       DO 300 I = 1,3
	       DO 300 J = 1,3
	       DO 300 K = 1,3
	       DO 300 L = 1,3
	          SJMAT(I,J,K,L) = 0.0
	      	AUXMAT1(I,J,K,L) = 0.0
	       DO 298  N = 1,3		     
	          AUXMAT1(I,J,K,L) =  AUXMAT1(I,J,K,L) +
     &                 FTTAU(I,N)*AUXMAT2(N,J,K,L)
298          CONTINUE
                SJMAT(I,J,K,L) = SJMAT(I,J,K,L) 
     &            + RTTAU(I,K)*FSTART(L,J)
     &            - RTTAU(I,K)*AUX2(L,J)
     &            - AUXMAT1(I,J,K,L)
300          CONTINUE	
C
C            STEP J13. CALCULATE THE JACOBIAN MATRIX FOR EACH CRYSTAL
C
             CALL M3INV(FSTAU,FSTAUINV)
	       CALL MDET(FSTAU,DETFSTAU)
C	     
	       DO  400 I = 1,3
	       DO  400 J = 1,3
	       DO  400 K = 1,3
	       DO  400 L = 1,3
	           WJMAT(I,J,K,L) = 0.0
		         SUM = 0.0
		     DO 350 IP = 1,3
		     DO 350 IQ = 1,3
		         SUM = SUM + SJMAT(IP,IQ,K,L)*FSTAUINV(IQ,IP)
350          CONTINUE
             DO 375 M = 1,3
		     DO 375 N = 1,3
		              WJMAT(I,J,K,L) = WJMAT(I,J,K,L) +
     &                SJMAT(I,M,K,L)*TSTARTAU(M,N)*FSTAU(J,N) +
     &                FSTAU(I,M)*QJMAT(M,N,K,L)*FSTAU(J,N) +
     &                FSTAU(I,M)*TSTARTAU(M,N)*SJMAT(J,N,K,L) -
     &                FSTAU(I,M)*TSTARTAU(M,N)*FSTAU(J,N)*SUM
375          CONTINUE
                WJMAT(I,J,K,L) = WJMAT(I,J,K,L)/DETFSTAU
390          CONTINUE
400          CONTINUE
C
C            STEP J14. CALCULATE THE REDUCED JACOBIAN MATRIX FOR EACH CRYSTAL
C
             CALL EMATR(WJMAT,WJMATR)
C
C          STEP J15. ADD THE COMPONENTS OF WJMATR FOR EACH CRYSTAL
C                    TO CALCULATE THE GLOBAL AVERAGED JACOBIAN
C
      	     DO 476 I = 1,6
	       DO 476 J = 1,6
	          WJMATGR(I,J) =  WJMATGR(I,J) + WJMATR(I,J)
476          CONTINUE	     		
C
       ITERK  = ITERKMAX

       RETURN
       END 
       
       
C*************************************************************************
	SUBROUTINE SOLVCRYS(GAMMAT,TSTV,ST,TSTRV,PMAT,CMAT,DTIME,TSTAUV,
	1           STAU,DGAMMA,DGMAX,RSS,DOMPSC,IITERK,IITERL,IITERERR,
     2		     RHO,RHO_NEW,ssd0,ssd_new,gnd0,gnd_new,NOEL)
      use modname
c	SUBROUTINE SOLVCRYS(GAMMAT,TSTV,ST,TSTRV,PMAT,CMAT,DTIME,TSTAUV,
c	1           STAU,DGAMMA,DGMAX,RSS,DOMPSC,IITERK,IITERL,IITERERR)

C	THIS SUBROUTINE SOLVES FOR  (TSTARTAU,STAU)
C       HERE, WE ARE SETTING
C
C       TSTARTV   AS TSTV
C       TSTARTRV  AS TSTRV
C       TSTARTAUV AS TSTAUV
C --------------------------------------------------------------------

 	IMPLICIT REAL*8(A-H,O-Z)
	
!	Include 'D:\temp-wagoner102\COMMON'
       
      PARAMETER  (MAXIT1=100,MAXIT2=15,NP=6,N=6)	  
      DIMENSION  ALPHAM(NP,NP),BETA(NP),INDX(NP)
	DIMENSION  TSTV(6),ST(MAXSLIP),TSTRV(6)
	DIMENSION  PMAT(MAXSLIP,6)
	DIMENSION  CMAT(MAXSLIP,6)
	DIMENSION  TSTAUV(6),STAU(MAXSLIP),STRESS(6)
	DIMENSION  RSS(MAXSLIP),DGAMMA(MAXSLIP),DGAMDTAU(MAXSLIP)
	DIMENSION  HAB(MAXSLIP,MAXSLIP)
	DIMENSION  RESOLD(MAXSLIP),RESNEW(MAXSLIP)
	DIMENSION  TEMPSTRESS(6)
	DIMENSION  GAMMAT(MAXSLIP)
	dimension  RHO(MAXSLIP),RHO_NEW(MAXSLIP),drho(maxslip)
	dimension  rho_old(maxslip),ssd_old(maxslip),gnd_old(maxslip)
	dimension  ssd0(maxslip),ssd_new(maxslip),gnd0(maxslip),
	1gnd_new(maxslip)
	
	IITERK    = 0
	IITERL    = 0
	IITERLNEW = 0
      TOL      = S0*1.E-4
	DGMAXTOL = 50000.0
C
C	INITIAL ESTIMATES FOR THE ROOTS:
C
	DO 10 I=1,6
	   STRESS(I) = TSTV(I)
10	CONTINUE

	DO 15 ISLIP=1,NSLIP
	   RESOLD(ISLIP) = ST(ISLIP)
15	CONTINUE

cnew

	do islip=1,nslip
		rho_old(islip)= RHO(islip)
		ssd_old(islip)= ssd0(islip)
		gnd_old(islip)= gnd0(islip)
	enddo
cnew

        
999   CONTINUE
      IITERERR = 0
C
C		START THE LOOP FOR SOLVING FOR THE COMPONENTS OF TSTAUV
C
	DO 500  K=1,MAXIT1
C
C	CALCULATE:
C     1. THE RESOLVED SHEAR STRESS ON EACH
C	  SLIP SYSTEM CORRESPONDING TO THE ESTIMATE OF TSTAUV.
C     2. CALCULATE THE PLASTIC SHEAR STRAIN INCREMENT FOR
C	  EACH SLIP SYSTEM.
C     3. CALCULATE THE MAXIMUM SHEAR STRAIN INCREMENT 
C     4. CALCULATE THE PLASTIC WORK INCREMENT
C
	 DGMAX   = 0.0
       DOMPSC = 0.0	
	 DO 30 ISLIP = 1,NSLIP
	     TAUALPHA = 0.0
	 DO 20 J=1,6
	     TAUALPHA = TAUALPHA + STRESS(J)*PMAT(ISLIP,J)
20	 CONTINUE

                RSS(ISLIP) = TAUALPHA
                RESALPHA   = RESOLD(ISLIP)

	CALL GAMMADOT(TAUALPHA,RESALPHA,GDOTALPHA,DGDTAU,NOEL,ISLIP
	1,temptau,kinc)	

c		rss(islip)=taualpha
		rss(islip)=temptau
					

						
        DGAMMA(ISLIP) = DTIME*GDOTALPHA
        DGAMDTAU(ISLIP) = DTIME*DGDTAU
		DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
		TEMP =  DABS ( RSS(ISLIP)*DGAMMA(ISLIP)) 
	
        DOMPSC = DOMPSC + TEMP
30	  CONTINUE
C
C		Calculate BETA_i --- the negative of the values of
C		the functions G_i(T_k)
C
	   DO 40 I =1,6
	      BETA(I) = 0.0
	      DO 35 ISLIP=1,NSLIP
	         BETA(I) = BETA(I) + DGAMMA(ISLIP)*CMAT(ISLIP,I)
35		 CONTINUE
	         BETA(I) =   TSTRV(I) - BETA(I) - STRESS(I)
40	   CONTINUE
C
C	       Calculate the Jacobian for the Newton corrections
C
	   DO 70 I =1,6
	   DO 70 J =1,6
               IF (I .EQ. J) THEN      
	            ALPHAM(I,J)      = 1.0D0
	         ELSE
	            ALPHAM(I,J)      = 0.0D0
	         END IF
	         DO 60 ISLIP=1,NSLIP
                     IF(RSS(ISLIP).EQ.0.0)GO TO 60
	             ALPHAM(I,J)= ALPHAM(I,J) + 
     &                     DGAMDTAU(ISLIP)*CMAT(ISLIP,I)*PMAT(ISLIP,J)
60	   CONTINUE
70	   CONTINUE    
C
C	         Calculate the Newton corrections
C
          CALL LUDCMP(ALPHAM,N,NP,INDX,D)
          CALL LUBKSB(ALPHAM,N,NP,INDX,BETA)
C
C	    BETA  NOW CONTAINS THE NEWTON CORRECTIONS.
C		  CHECK TO SEE IF THE NEWTON CORRECTIONS ARE
C	    WITHIN SPECIFIED TOLERANCES. IF YES, THEN
C		  PROCESS HAS CONVERGED. UPDATE AND GO TO LOOP FOR
C         ITERATING ON THE DEFORMATION RESISTANCES.
C		  IF NOT, THEN CONTINUE ITERATING.
          IF (  ( DABS( BETA(1) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(2) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(3) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(4) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(5) ) .LE. TOL ) .AND.
     +	         ( DABS( BETA(6) ) .LE. TOL)   )  THEN  
C
C		     ITERATIVE PROCEDURE HAS CONVERGED. UPDATE 
C
	      TSTAUV(1) = STRESS(1)
	      TSTAUV(2) = STRESS(2)
	      TSTAUV(3) = STRESS(3)
	      TSTAUV(4) = STRESS(4)
	      TSTAUV(5) = STRESS(5)
	      TSTAUV(6) = STRESS(6)	      
            IITERK    =  K
              GO TO 550
           END IF
C
C		TEMPORARILY ADD THE CORRECTIONS TO THE VARIABLES
C
	     TEMPSTRESS(1) = STRESS(1) + BETA(1)
	     TEMPSTRESS(2) = STRESS(2) + BETA(2)
	     TEMPSTRESS(3) = STRESS(3) + BETA(3)
	     TEMPSTRESS(4) = STRESS(4) + BETA(4)
	     TEMPSTRESS(5) = STRESS(5) + BETA(5)
	     TEMPSTRESS(6) = STRESS(6) + BETA(6)
C
C               CHECK IF THESE STRESS COMPONENTS GIVE RISE TO LARGE DELTA GAMMAS
C               IF THEY DO THEN REDUCE THE N-R CORRECTIONS BY A FACTOR OF 4
C               AND CONTINUE UNTIL THE DELTA GAMMAS ARE REASONABLE
C
           ICOR = 0
85         CONTINUE
           IF(ICOR .GT. 5) THEN 
C             THE METHOD TO CONSTRAIN THE NEWTON-RAPHSON METHOD 
C             IS NOT WORKING! GET OUT OF THE LOOP! 
C             RETURN TO CUT BACK TIME-STEP AND RECALCULATE
C 
          if (imeg.eq.1) WRITE(UR7,*) 'TOO MANY ITERATIONS IN CONSTRAINT'

	     GO TO 500
           END IF
c---------------------------------
	  DGMAX = 0.0
	  DO 100 ISLIP = 1,NSLIP
	        TAUALPHA = 0.0
		DO 90 J=1,6
		     TAUALPHA = TAUALPHA + TEMPSTRESS(J)*PMAT(ISLIP,J)
90		CONTINUE
	

                RSS(ISLIP) = TAUALPHA
                RESALPHA   = RESOLD(ISLIP)

         	
	CALL GAMMADOT(TAUALPHA,RESALPHA,GDOTALPHA,DGDTAU,NOEL,ISLIP
	1,temptau,kinc)	
      
c		  RSS(ISLIP) = TAUALPHA
		  RSS(ISLIP) = Temptau
		  

			  
		  DGAMMA(ISLIP) = DTIME*GDOTALPHA
          DGAMDTAU(ISLIP) = DTIME*DGDTAU
	    DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
100	   CONTINUE
C
           IF (DGMAX .GT. DGMAXTOL) THEN
	      DO 102 I = 1,6
	         BETA(I) = 0.25 * BETA(I)
		       TEMPSTRESS(I) = STRESS(I) + BETA(I)
102           CONTINUE
              ICOR = ICOR +1
	
		    GO TO 85
	   ELSE

C
C	       ACCEPT THE CORRECTIONS TO THE STRESS AND CONTINUE ITERATING
C
              DO 104 I = 1,6
	         STRESS(I) = TEMPSTRESS(I) 
104           CONTINUE
           END IF
C	    	     
500	CONTINUE
C			ITERATIVE PROCEDURE FOR THE STRESS 
C			HAS NOT CONVERGED IN MAXIT1 ITERATIONS.
C			RETURN TO CUT BACK TIME-STEP, AND RECALCULATE.
C
      IITERK = MAXIT1
	IITERERR = 1

	RETURN
C
C
550   CONTINUE

C		ITERATIVE PROCEDURE FOR TSTAUV HAS CONVERGED
C		START ITERATING FOR THE SLIP SYSTEM DEFORMATION 
C		RESISTANCES

C
C		ESTIMATE THE THE SLIP SYSTEM DEFORMATION RESISTANCES
C		AT THE END OF THE STEP

c------------------------------------------------------------------
	call hardening3(NOEL,rho_old,dgamma,dtime,hab,drho)
c------------------------------------------------------------------

c	INPUT: old dens,shear increment, OUTPUT: hab, drho	


C		ESTIMATE THE SLIP SYSTEM DEFORMATION RESISTANCES AT THE
C		END OF THE INCREMENT
C	
		! new hardeness based on disloction density
		! April: decompose SSD and GND
c-----------------------------------------------------------------------
      
	do islip=1,nslip
	SSD_NEW(ISLIP)=SSD_OLD(ISLIP)+DRHO(ISLIP)
	gnd_new(islip)=gnd_old(islip)
	enddo
	     
C---------------------------------------------------------------------------
	do islip=1,nslip
	rho_new(islip)=ssd_new(islip)+ABS(gnd_new(islip))
	enddo 
c--------------------------------------------------------------------		
		
	do islip=1,nslip
	tempres=0.
		do jslip=1,nslip
		tempres=tempres+hab(islip,jslip)*rho_new(jslip)
		enddo
		RESNEW(ISLIP)=0.4*SHM*BS*DSQRT(TEMPRES)	
c      RESNEW(ISLIP)=10
	enddo

	ERRRES = 0.0
	DO 650 ISLIP = 1,NSLIP
         ERRRES = DMAX1(DABS(RESNEW(ISLIP)-RESOLD(ISLIP)),ERRRES)
650   CONTINUE

      TOLS = S0*1.D-3
      IF (ERRRES .LT. TOLS) THEN
C		ITERATIVE PROCEDURE HAS CONVERGED
C		UPDATE STAU AND RETURN
C
         DO 660 ISLIP = 1,NSLIP
            STAU(ISLIP) = RESNEW(ISLIP)
660      CONTINUE
         RETURN
      ELSE
         DO 670 ISLIP = 1,NSLIP
            RESOLD(ISLIP) = RESNEW(ISLIP)
670   CONTINUE
         IITERLNEW = IITERLNEW + 1
         IITERL = MAX0(IITERL ,IITERLNEW)
	IF (ITERL .GT. MAXIT2) THEN
C		ITERATIVE PROCEDURE FOR THE DEFORMATION RESISTANCE HAS NOT
C		CONVERGED IN MAXIT2 ITERATIONS. RETURN TO CUT BACK THE TIME
C		STEP, AND RECALCULATE.
	   IITERL = MAXIT2
	   IITERERR = 1
	if (imeg.eq.1) WRITE (UR7,*) 'MAXIT2 EXCEEDED....REPEATING TIME STEP'
	RETURN
	ENDIF
      GO TO 999
      END IF

      END

c--------------------------------------------------------------------------

	SUBROUTINE GAMMADOT(TAUALPHA,RESALPHA,GDOTALPHA,
	1 DGDTAU,NOEL,ISLIP,temptau,kinc)

      use modname
	IMPLICIT REAL*8(A-H,O-Z) 
C	real(8) bftemp,TEMPTAU,TAUOBS
C	integer noel,islip,i,kinc
C      DIMENSION htrans1(maxnoel,maxslip)

!      Include 'D:\temp-wagoner102\COMMON'		
      
c--------------------------------------------------------------------------

	IF (taualpha.eq.0) then
	bftemp=0
	ELSE
	bftemp=abs(bftt(noel,islip))
	
	ENDIF
	
	
c	write(UR10,*)resalpha

c-------Transmission criteria----------------------------------------------	
		
	IF (GTYP(NOEL,ISLIP).EQ.2.and.taualpha.ne.0) THEN
c-twq
	tauobs=taustar*(1-htrans(noel,islip)) 

	bftemp=tauobs+bftemp

		
	ELSE
		
	ENDIF

c--------Back Stress--------------------------------------------------------	

	IF (taualpha.ge.0) THEN
	temptau=taualpha-bftemp
	
		IF (temptau.LE.0) then
		TEMPTAU=0
	  GDOTALPHA = 0.0
		DGDTAU   = 0.0
		goto 503		
		ENDIF

	ELSE							!TAUALPHA<0
	
	temptau=taualpha+bftemp

		IF (temptau.GE.0) THEN	
		TEMPTAU=0
		GDOTALPHA = 0.0
		DGDTAU   = 0.0
		goto 503
		ENDIF

	endif	

c----------------------------------------------------------------------------
	
      IF(DABS(temptau) .LE. S0*1.E-3) THEN  
      GDOTALPHA = 0.0
	DGDTAU   = 0.0
	ELSE
	   IF (temptau.GE.0.D0) THEN
	   GDOTALPHA = GDOT0*DABS(temptau/RESALPHA)**(1./AM)
	   ELSE
		 GDOTALPHA = -GDOT0*DABS(temptau/RESALPHA)**(1./AM)
	   ENDIF
	   DGDTAU = GDOTALPHA/(AM*temptau)
      ENDIF

503	continue


      RETURN
      END



c--------------------------------------------------------------------------	
	
	SUBROUTINE HARDENING3(NOEL,RHO,DGAMMA,DTIME,HAB,DRHO)
            use modname	
	IMPLICIT REAL*8(A-H,O-Z)

!	Include 'D:\temp-wagoner102\COMMON'

	DIMENSION RHO(MAXSLIP),DGAMMA(MAXSLIP),HAB(MAXSLIP,MAXSLIP)
	DIMENSION DRHO(MAXSLIP)
C
C	CALCULATE THE HARDENING MATRIX HAB
C
	DO 10 ISLIP = 1,NSLIP

	     IF( (DABS(DGAMMA(ISLIP)) / DTIME) .LT. 1.E-10) THEN
	        
			DRHO(islip)=0.

	     ELSE
			temp=0.
			do jslip=1,nslip
				if(jslip.ne.islip) temp=temp+RHO(jslip)
				if(nslip.eq.1) temp=RHO(jslip)
			enddo
			
			if(temp.eq.0.) then
				
				temp1=0.
				do jslip=1,nslip
					if(jslip.ne.islip) temp1=temp1+rh0
				enddo

				temp=temp1

			endif
			
			temp0=clength(noel,islip)
			
			if(itype.eq.0.and.dk.eq.0.) then
				drho(islip)=0.
				goto 345
			endif

			if(itype.eq.0) then           ! itype=0
				pathL=dk*temp**(-1./2.)
			else if(itype.eq.1) then     ! itype 1
				pathL=clength(noel,islip)

			endif
			
			term1=dabs(DGAMMA(islip))/bs/pathL*decmp_mech
c			TERM1=0
			
			term2=dabs(DGAMMA(islip))/bs*
	1		(yc*rho(islip))*decmp_mech
c	      TERM2=0

			DRHO(islip)=term1-term2
c            DRHO(ISLIP)=0
			if(term2.gt.term1) then
				if (icheck.eq.1) write(UR8,*) 'drho<0'
				drho(islip)=0.
			endif

345			continue	

	     END IF
10	CONTINUE
      
C      WRITE(UR10,*)RH0,TEMP,PATHL

C
C    CALCULATE THE LATENT HARDENING MATRIX QLAT
C

100	continue

	DO ISLIP = 1,NSLIP
		DO JSLIP = 1,NSLIP
			if(islip.eq.jslip) then
				HAB(ISLIP,JSLIP)=ass
			else
				hab(islip,jslip)=asu
			endif
		enddo
	enddo


      RETURN 
      END

C***********************************************************************
	SUBROUTINE ROTMAT(TH,PHI,OM,Q)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION Q(3,3)
C
C		INITIALIZE
C
	CALL ONEM(Q)
        STH = DSIN(TH)
        CTH = DCOS(TH)
        SPH = DSIN(PHI)
        CPH = DCOS(PHI)
	  SOM = DSIN(OM)
        COM = DCOS(OM)
C
        Q(1,1) = CPH*COM-SPH*SOM*CTH
        Q(1,2) = SPH*COM+CPH*SOM*CTH
        Q(1,3) = SOM*STH
C
        Q(2,1) = -CPH*SOM-SPH*COM*CTH
        Q(2,2) = -SPH*SOM+CPH*COM*CTH
        Q(2,3) =  COM*STH
C
        Q(3,1) = SPH*STH
        Q(3,2) = -CPH*STH
        Q(3,3) = CTH
C
	RETURN
	END	        
C ****************************************************************
	SUBROUTINE ELASTOP(E,Q,QT,T)
            use modname
C
C		THIS SUBROUTINE CONVERTS A SYMMETRIC STRAIN
C		MATRIX E GIVEN IN GLOBAL COORDINATES TO
C		A MATRIX EC IN CRYSTAL COORDINATES. THEN CALCULATES
C		TC=L[EC] WITH L THE ANISOTROPIC ELASTIC STIFFNESS
C		MATRIX IN CRYSTAL COORDINATES FOR A CUBIC CRYSTAL,
C		AND FINALLY CONVERTS THE STRESS MATYRIX TC BACK TO
C		A MATRIX T IN GLOBAL COORDINATES
C
	IMPLICIT REAL*8(A-H,O-Z)
	
!      Include 'D:\temp-wagoner102\COMMON'

     

	DIMENSION E(3,3),Q(3,3),QT(3,3),EC(3,3),TC(3,3),T(3,3),AUX1(3,3)

	CALL MPROD(QT,E,AUX1)
	CALL MPROD(AUX1,Q,EC)
        TC(1,1)= C11*EC(1,1)+C12*EC(2,2)+C12*EC(3,3)
        TC(2,2)= C12*EC(1,1)+C11*EC(2,2)+C12*EC(3,3)        
        TC(3,3)= C12*EC(1,1)+C12*EC(2,2)+C11*EC(3,3)       
        TC(1,2)=2.*C44*EC(1,2)
        TC(1,3)=2.*C44*EC(1,3)
        TC(2,3)=2.*C44*EC(2,3)
        TC(2,1)=TC(1,2)
        TC(3,1)=TC(1,3)
        TC(3,2)=TC(2,3)
	CALL MPROD(Q,TC,AUX1)
	CALL MPROD(AUX1,QT,T)

	RETURN
	END
C*********************************************************************
      SUBROUTINE ELASMAT(Q,ELMAT)
            use modname
C      
C          THIS SUBROUTINE CALCULATES THE ELASTIC STIFFNESS MATRIX 
C          AS A 3 X 3 X 3 X 3 MATRIX IN GLOBAL COORDINATES
C
	IMPLICIT REAL*8(A-H,O-Z)
	
!      Include 'D:\temp-wagoner102\COMMON'

      
      DIMENSION Q(3,3),ELMAT(3,3,3,3)
      DIMENSION ELMATC(3,3,3,3)
      DIMENSION AUXMAT1(3,3,3,3),AUXMAT2(3,3,3,3)
      
      CALL ZEROMOD(ELMATC)
      ELMATC(1,1,1,1) = C11
      ELMATC(1,1,2,2) = C12
      ELMATC(1,1,3,3) = C12
      ELMATC(2,2,1,1) = C12
      ELMATC(2,2,2,2) = C11
      ELMATC(2,2,3,3) = C12      
      ELMATC(3,3,1,1) = C12  
      ELMATC(3,3,2,2) = C12      
      ELMATC(3,3,3,3) = C11
      ELMATC(1,2,1,2) = C44
      ELMATC(1,2,2,1) = C44      
      ELMATC(2,1,1,2) = C44  
      ELMATC(2,1,2,1) = C44
      ELMATC(1,3,1,3) = C44
      ELMATC(1,3,3,1) = C44
      ELMATC(3,1,1,3) = C44    
      ELMATC(3,1,3,1) = C44   
      ELMATC(2,3,2,3) = C44  
      ELMATC(2,3,3,2) = C44 
      ELMATC(3,2,2,3) = C44  
      ELMATC(3,2,3,2) = C44 
  
       DO 20 I=   1,3
       DO 20 N =  1,3
       DO 20 IO = 1,3
       DO 20 IP = 1,3
             AUXMAT1(I,N,IO,IP) = 0.0
	     DO 15  M = 1,3
	        AUXMAT1(I,N,IO,IP) = AUXMAT1(I,N,IO,IP) +
     &                              ELMATC(M,N,IO,IP)*Q(I,M)
15           CONTINUE
20     CONTINUE

       DO 30 I=   1,3
       DO 30 J =  1,3
       DO 30 IO = 1,3
       DO 30 IP = 1,3
             AUXMAT2(I,J,IO,IP) = 0.0
	     DO 25  N = 1,3
	        AUXMAT2(I,J,IO,IP) = AUXMAT2(I,J,IO,IP) +
     &                              AUXMAT1(I,N,IO,IP)*Q(J,N)
25           CONTINUE
30     CONTINUE

       DO 40 I  = 1,3
       DO 40 J  = 1,3
       DO 40 K  = 1,3
       DO 40 IP = 1,3
             AUXMAT1(I,J,K,IP) = 0.0
	     DO 35  IO = 1,3
	        AUXMAT1(I,J,K,IP) = AUXMAT1(I,J,K,IP) +
     &                              AUXMAT2(I,J,IO,IP)*Q(K,IO)
35           CONTINUE
40     CONTINUE

       DO 50 I  = 1,3
       DO 50 J  = 1,3
       DO 50 K  = 1,3
       DO 50 L  = 1,3
             ELMAT(I,J,K,L) = 0.0
	     DO 45  IP = 1,3
	        ELMAT(I,J,K,L) = ELMAT(I,J,K,L) +
     &                              AUXMAT1(I,J,K,IP)*Q(L,IP)
45           CONTINUE
50     CONTINUE
       
       RETURN
       END	
C**********************************************************************
      SUBROUTINE EMATR(ELMAT,ELMATR)
C      
C          THIS SUBROUTINE REDUCES AN ELASTIC STIFFNESS TYPE MATRIX 
C          3 X 3 X 3 X 3 MATRIX TO A 6 X 6 MATRIX
C          
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMAT(3,3,3,3), ELMATR(6,6)
      
        DO 10 I = 1,6
	DO 10 J = 1,6
	    ELMATR(I,J) = 0.0
10      CONTINUE
      
 
        ELMATR(1,1)  = ELMAT(1,1,1,1)
	ELMATR(1,2)  = ELMAT(1,1,2,2)
	ELMATR(1,3)  = ELMAT(1,1,3,3) 
	ELMATR(1,4)  = 0.5 * ( ELMAT(1,1,1,2) + ELMAT(1,1,2,1) )
	ELMATR(1,5)  = 0.5 * ( ELMAT(1,1,1,3) + ELMAT(1,1,3,1) )
	ELMATR(1,6)  = 0.5 * ( ELMAT(1,1,2,3) + ELMAT(1,1,3,2) )


        ELMATR(2,1)  = ELMAT(2,2,1,1)
	ELMATR(2,2)  = ELMAT(2,2,2,2)
	ELMATR(2,3)  = ELMAT(2,2,3,3) 
	ELMATR(2,4)  = 0.5 * ( ELMAT(2,2,1,2) + ELMAT(2,2,2,1) )
	ELMATR(2,5)  = 0.5 * ( ELMAT(2,2,1,3) + ELMAT(2,2,3,1) )
	ELMATR(2,6)  = 0.5 * ( ELMAT(2,2,2,3) + ELMAT(2,2,3,2) )
       
        ELMATR(3,1)  = ELMAT(3,3,1,1)
	ELMATR(3,2)  = ELMAT(3,3,2,2)
	ELMATR(3,3)  = ELMAT(3,3,3,3) 
	ELMATR(3,4)  = 0.5 * ( ELMAT(3,3,1,2) + ELMAT(3,3,2,1) )
	ELMATR(3,5)  = 0.5 * ( ELMAT(3,3,1,3) + ELMAT(3,3,3,1) )
	ELMATR(3,6)  = 0.5 * ( ELMAT(3,3,2,3) + ELMAT(3,3,3,2) )   
	
        ELMATR(4,1)  = ELMAT(1,2,1,1)
	ELMATR(4,2)  = ELMAT(1,2,2,2)
	ELMATR(4,3)  = ELMAT(1,2,3,3) 
	ELMATR(4,4)  = 0.5 * ( ELMAT(1,2,1,2) + ELMAT(1,2,2,1) )
	ELMATR(4,5)  = 0.5 * ( ELMAT(1,2,1,3) + ELMAT(1,2,3,1) )
	ELMATR(4,6)  = 0.5 * ( ELMAT(1,2,2,3) + ELMAT(1,2,3,2) ) 	
	
	ELMATR(5,1)  = ELMAT(1,3,1,1)
	ELMATR(5,2)  = ELMAT(1,3,2,2)
	ELMATR(5,3)  = ELMAT(1,3,3,3) 
	ELMATR(5,4)  = 0.5 * ( ELMAT(1,3,1,2) + ELMAT(1,3,2,1) )
	ELMATR(5,5)  = 0.5 * ( ELMAT(1,3,1,3) + ELMAT(1,3,3,1) )
	ELMATR(5,6)  = 0.5 * ( ELMAT(1,3,2,3) + ELMAT(1,3,3,2) ) 

	ELMATR(6,1)  = ELMAT(2,3,1,1)
	ELMATR(6,2)  = ELMAT(2,3,2,2)
	ELMATR(6,3)  = ELMAT(2,3,3,3) 
	ELMATR(6,4)  = 0.5 * ( ELMAT(2,3,1,2) + ELMAT(2,3,2,1) )
	ELMATR(6,5)  = 0.5 * ( ELMAT(2,3,1,3) + ELMAT(2,3,3,1) )
	ELMATR(6,6)  = 0.5 * ( ELMAT(2,3,2,3) + ELMAT(2,3,3,2) ) 
	
       
        RETURN
        END
C******************************************************************************
      SUBROUTINE EMATREC(ELMATR,ELMAT)
C      
C          THIS SUBROUTINE RECONSTRUCTS A STIFFNESS MATRIX ELMATREC
C          AS A 3 X 3 X 3 X 3 MATRIX FROM ITS 
C          REDUCED 6 X 6 FORM ELMATR
C
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMATR(6,6),ELMAT(3,3,3,3)
      
      DO 10 I = 1,3
      DO 10 J = 1,3
      DO 10 K = 1,3
      DO 10 L = 1,3
           ELMAT(I,J,K,L) = 0.0
10    CONTINUE
      

C
C  RECONSTRUCT THE MATRIX ELMAT FROM ITS REDUCED FORM ELMATR
C
	  
        ELMAT(1,1,1,1)  = ELMATR(1,1)
	ELMAT(1,1,2,2)  = ELMATR(1,2)
	ELMAT(1,1,3,3)  = ELMATR(1,3) 
	ELMAT(1,1,1,2)  = ELMATR(1,4)
	ELMAT(1,1,2,1)  = ELMATR(1,4)
	ELMAT(1,1,1,3)  = ELMATR(1,5)
	ELMAT(1,1,3,1)  = ELMATR(1,5)
	ELMAT(1,1,2,3)  = ELMATR(1,6)
	ELMAT(1,1,3,2)  = ELMATR(1,6)
	

        ELMAT(2,2,1,1)  = ELMATR(2,1)
	ELMAT(2,2,2,2)  = ELMATR(2,2)
	ELMAT(2,2,3,3)  = ELMATR(2,3) 
	ELMAT(2,2,1,2)  = ELMATR(2,4)
	ELMAT(2,2,2,1)  = ELMATR(2,4)
	ELMAT(2,2,1,3)  = ELMATR(2,5)
	ELMAT(2,2,3,1)  = ELMATR(2,5)
	ELMAT(2,2,2,3)  = ELMATR(2,6)
	ELMAT(2,2,3,2)  = ELMATR(2,6)
	
	ELMAT(3,3,1,1)  = ELMATR(3,1)
	ELMAT(3,3,2,2)  = ELMATR(3,2)
	ELMAT(3,3,3,3)  = ELMATR(3,3) 
	ELMAT(3,3,1,2)  = ELMATR(3,4)
	ELMAT(3,3,2,1)  = ELMATR(3,4)
	ELMAT(3,3,1,3)  = ELMATR(3,5)
	ELMAT(3,3,3,1)  = ELMATR(3,5)
	ELMAT(3,3,2,3)  = ELMATR(3,6)
	ELMAT(3,3,3,2)  = ELMATR(3,6)
	
        ELMAT(1,2,1,1)  =  ELMATR(4,1)
        ELMAT(2,1,1,1)  =  ELMATR(4,1)
        ELMAT(1,2,2,2)  =  ELMATR(4,2)
        ELMAT(2,1,2,2)  =  ELMATR(4,2)	
	ELMAT(1,2,3,3)  =  ELMATR(4,3)
        ELMAT(2,1,3,3)  =  ELMATR(4,3)
	ELMAT(1,2,1,2)  =  ELMATR(4,4) 
	ELMAT(2,1,1,2)  =  ELMATR(4,4) 
        ELMAT(1,2,2,1)  =  ELMATR(4,4)
	ELMAT(2,1,2,1)  =  ELMATR(4,4)
	ELMAT(1,2,1,3)  =  ELMATR(4,5) 
	ELMAT(2,1,1,3)  =  ELMATR(4,5) 
        ELMAT(1,2,3,1)  =  ELMATR(4,5)
	ELMAT(2,1,3,1)  =  ELMATR(4,5)
	ELMAT(1,2,2,3)  =  ELMATR(4,6) 
	ELMAT(2,1,2,3)  =  ELMATR(4,6) 
        ELMAT(1,2,3,2)  =  ELMATR(4,6)
	ELMAT(2,1,3,2)  =  ELMATR(4,6)
	
	
	ELMAT(1,3,1,1)  =  ELMATR(5,1)
        ELMAT(3,1,1,1)  =  ELMATR(5,1)	
        ELMAT(1,3,2,2)  =  ELMATR(5,2)
        ELMAT(3,1,2,2)  =  ELMATR(5,2)	
	ELMAT(1,3,3,3)  =  ELMATR(5,3)
        ELMAT(3,1,3,3)  =  ELMATR(5,3)
	ELMAT(1,3,1,2)  =  ELMATR(5,4) 
	ELMAT(3,1,1,2)  =  ELMATR(5,4) 
        ELMAT(1,3,2,1)  =  ELMATR(5,4)
	ELMAT(3,1,2,1)  =  ELMATR(5,4)
	ELMAT(1,3,1,3)  =  ELMATR(5,5) 
	ELMAT(3,1,1,3)  =  ELMATR(5,5) 
        ELMAT(1,3,3,1)  =  ELMATR(5,5)
	ELMAT(3,1,3,1)  =  ELMATR(5,5)
	ELMAT(1,3,2,3)  =  ELMATR(5,6) 
	ELMAT(3,1,2,3)  =  ELMATR(5,6) 
        ELMAT(1,3,3,2)  =  ELMATR(5,6)
	ELMAT(3,1,3,2)  =  ELMATR(5,6)
	
	ELMAT(2,3,1,1)  =  ELMATR(6,1)
        ELMAT(3,2,1,1)  =  ELMATR(6,1)	
        ELMAT(2,3,2,2)  =  ELMATR(6,2)
        ELMAT(3,2,2,2)  =  ELMATR(6,2)	
	ELMAT(2,3,3,3)  =  ELMATR(6,3)
        ELMAT(3,2,3,3)  =  ELMATR(6,3)
	ELMAT(2,3,1,2)  =  ELMATR(6,4) 
	ELMAT(3,2,1,2)  =  ELMATR(6,4) 
        ELMAT(2,3,2,1)  =  ELMATR(6,4)
	ELMAT(3,2,2,1)  =  ELMATR(6,4)
	ELMAT(2,3,1,3)  =  ELMATR(6,5) 
	ELMAT(3,2,1,3)  =  ELMATR(6,5) 
        ELMAT(2,3,3,1)  =  ELMATR(6,5)
	ELMAT(3,2,3,1)  =  ELMATR(6,5)
	ELMAT(2,3,2,3)  =  ELMATR(6,6) 
	ELMAT(3,2,2,3)  =  ELMATR(6,6) 
        ELMAT(2,3,3,2)  =  ELMATR(6,6)
	ELMAT(3,2,3,2)  =  ELMATR(6,6)

       RETURN 
       END  
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a row-wise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the nuber of
C	row interchanges was even or odd, respectively. This routine
C	is used in combination with LUBKSB to solve linear equations 
C	or invert a matrix.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	
	PARAMETER (NMAX=1000,TINY=1.0E-20)
	DIMENSION A(NP,NP),INDX(N),VV(N)
	D=1.
	DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) then
	  write(*,*) 'Singular matrix.'
		stop
		endif
		
        VV(I)=1./AAMAX
12	CONTINUE
	DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16	CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19	CONTINUE
	IF(A(N,N).EQ.0.)A(N,N)=TINY
	RETURN
	END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C	Solves the set of N linear equations [A]{X} = {B}. 
C	Here [A] is input, not as the matrix [A], but as its LU 
C	decomposition, determined by the routine LUDCMP. INDX
C	is input as the permutation vector returned by LUDCMP. {B}
C	is input as the right-hand side vector {B}, and returns
C	with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and can be left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12     CONTINUE
       DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14     CONTINUE
       RETURN
       END
C ****************************************************************
      SUBROUTINE MATINV(A,N,NP,INDX,Y)

C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a rowwise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the number of
C	row interchanges was even or odd, respectively.
C
C	Once the LU decomposition is performed, this routine 
C	calculates the inverse of [A] by using subroutine LUBKSB.
C	Note that INDX is input as the permutation vector returned by 	
C	LUDCMP. {B} is input as the right-hand side vector {B}, and 	
C	returns with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and are left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
C	The inverse of [A] is calculated using as many unit vectors
C	{B} needed as right hand side vectors. The result is
C	returned as the matrix [Y].
C
	IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), Y(NP,NP), INDX(NP)
C
C	Set up the identity matrix  
C
	DO 12 I = 1,N
	    DO 11  J = 1,N
	        Y(I,J) = 0.0
11	    CONTINUE
          Y(I,I) = 1.
12	CONTINUE
C
C	Decompose the matrix just once
C
	CALL LUDCMP(A,N,NP,INDX,D)
C
C	Find the inverse by columns. It is necessary to recognize
C	that FORTRAN stores two dimensional matrices by column, so
C	so that Y(1,J) is the address of the Jth column of Y.
C
	DO 13 J=1,N
	    CALL LUBKSB(A,N,NP,INDX,Y(1,J))
13	CONTINUE
	RETURN
	END
C ********************************************************************
C
C	THE NEXT SUBROUTINE CALCULATES VARIOUS KINEMATICAL QUANTITIES 
C	ASSOCIATED WITH THE DEFORMATION GRADIENT
C
C *******************************************************************
	SUBROUTINE SKINEM(F,R,U)
C
C 	SUBROUTINE SKINEM(F,R,U,UEIGVAL,EIGVEC,E)
C
C		THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
C		[F]=[R][U] OF THE DEFORMATION GRADIENT [F] INTO
C		A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
C		THE EIGENVALUES AND EIGENVECTORS OF [U] AND
C		THE  LOGARITHMIC STRAIN [E] = LN [U]
C		ARE ALSO RETURNED.
C ------------------------------------------------------------------
C	VARIABLES
C
	REAL*8 F(3,3),DETF,FTRANS(3,3),
     +       C(3,3), OMEGA(3),
     +	 UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +	 U(3,3),E(3,3),UINV(3,3),R(3,3),TEMP(3,3)
         
C
C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C			   POLAR DECOMPOSITION IS DESIRED.
C	DETF		-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C			   OUTPUT.
C	U(3,3)	-- THE RIGHT STRETCH TENSOR; SYMMETRIC
C			   AND POSITIVE DEFINITE; OUTPUT.
C	UINV(3,3)	-- THE INVERSE OF [U].
C	C(3,3)	-- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
C			   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)	-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	UEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
C	EIGVECT(3,3)-- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
C			   OUTPUT.
C -----------------------------------------------------------------
C	COMPUTATION
C	

C
C		STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]
C
	CALL ONEM(R)
	CALL ONEM(U)
	CALL ONEM(UINV)
C
C		STORE THE ZERO MATRIX IN [E]
C
	CALL ZEROM(E)

C
C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C		IF NOT, THEN PRINT DIAGONOSTIC AND STOP.
C
       CALL MDET(F,DETF)

       IF (DETF .LE. 0.D0) THEN
          if (imeg.eq.1) WRITE(UR7,100)
          RETURN
       END IF
C
C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]
C
       CALL  MTRANS(F,FTRANS)
       CALL  MPROD(FTRANS,F,C)

C 
C		CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]
C

	 CALL SPECTRAL(C,OMEGA,EIGVEC)

C
C		CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]
C
	
	 UEIGVAL(1) = DSQRT(OMEGA(1))
	 UEIGVAL(2) = DSQRT(OMEGA(2))
	 UEIGVAL(3) = DSQRT(OMEGA(3))

	 U(1,1) = UEIGVAL(1)
	 U(2,2) = UEIGVAL(2)
	 U(3,3) = UEIGVAL(3)

	 E(1,1) = DLOG( UEIGVAL(1) )
	 E(2,2) = DLOG( UEIGVAL(2) )
	 E(3,3) = DLOG( UEIGVAL(3) )

C
C	CALCULATE THE COMPLETE TENSORS [U] AND [E]
C

	CALL MTRANS(EIGVEC,EIGVECT)
	CALL MPROD(EIGVEC,U,TEMP)
	CALL MPROD(TEMP,EIGVECT,U)
	CALL MPROD(EIGVEC,E,TEMP)
	CALL MPROD(TEMP,EIGVECT,E)

C
C		CALCULATE [UINV]
C

	CALL  M3INV(U,UINV)

C
C		CACULATE [R]

C
	CALL MPROD(F,UINV,R)
C -------------------------------------------------------------------
C	FORMATS
C
 100  FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +       ' IS NOT GREATER THAN 0',/,
     +       10X,'POLAR DECOMPOSITION NOT PERFORMED')
C -------------------------------------------------------------------
	RETURN
	END
C ****************************************************************
C
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C
C **********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)

C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C ---------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
              E(I,J) = A(I,J)
 1	    CONTINUE
 2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
	END
C *********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C
C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.
C
C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C ----------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C
C	 INITIALIZE [V] TO THE IDENTITY MATRIX
C

	DO 12 IP = 1,N	
	      DO 11 IQ = 1,N
	            V(IP,IQ) = 0.D0
 11          CONTINUE
             V(IP,IP) = 1.D0
 12	CONTINUE

C
C	INITIALIZE B AND D TO THE DIAGONAL OF [A], AND Z TO ZERO. THE
C	VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)
C

	DO 13 IP = 1,N
	    B(IP) = A(IP,IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
 13	CONTINUE

	NROT = 0

	DO 24 I = 1,50

C
C		SUM OFF-DIAGONAL ELEMENTS
C
	    SM = 0.D0
	    DO 15 IP = 1, N-1
	         DO 14 IQ = IP + 1, N

	             SM = SM + DABS ( A(IP,IQ ))

 14            CONTINUE
 15       CONTINUE
C
C		IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C		WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C		UNDERFLOW.
C
	      IF( SM .EQ. 0.D0) RETURN
C
C		IF( SM .LT. 1.0D-15) RETURN
C
C		IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C		|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, C	C		SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
C
          IF( I .LT. 4) THEN
               TRESH = 0.2D0*SM/N**2
          ELSE
               TRESH = 0.D0
          END IF

          DO 22 IP = 1, N-1
               DO 21 IQ = IP+1,N
                    G = 100.D0*DABS(A(IP,IQ))
C
C				AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C				OFF-DIAGONAL ELEMENT IS SMALL.
C

	IF( (I .GT. 4) .AND. ( DABS(D(IP))+G  .EQ.  DABS( D(IP)) )
     +               .AND. ( DABS(D(IQ))+G  .EQ.  DABS( D(IQ)) ) )THEN

                       A(IP,IQ) = 0.D0
                   ELSE IF( DABS(A(IP,IQ)) .GT. TRESH) THEN
                       H = D(IQ) - D(IP)
                       IF(DABS(H)+G .EQ. DABS(H)) THEN
C
C				T= 1./(2.*THETA), EQUATION(11.1.10)
C
			        T =A(IP,IQ)/H
			      ELSE
			          THETA = 0.5D0*H/A(IP,IQ)
			     T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
			          IF(THETA .LT. 0.D0) T = -T
			      END IF
			      C = 1.D0/DSQRT(1.D0 + T**2)
			      S = T*C
			      TAU = S/(1.D0 + C)
				H = T*A(IP,IQ)
				Z(IP) = Z(IP) - H
				Z(IQ) = Z(IQ) + H
				D(IP) = D(IP) - H
				D(IQ) = D(IQ) + H
				A(IP,IQ) = 0.D0
C
C				CASE OF ROTATIONS 1<= J < P
C				
				DO 16 J = 1, IP-1
				     G = A(J,IP)
				     H = A(J,IQ)
				     A(J,IP) = G - S*(H + G*TAU)
				     A(J,IQ) = H + S*(G - H*TAU)
 16				CONTINUE
C
C				CASE OF ROTATIONS P < J < Q
C
				DO 17 J = IP+1, IQ-1
				     G = A(IP,J)
				     H = A(J,IQ)
				     A(IP,J) = G - S*(H + G*TAU)
				     A(J,IQ) = H + S*(G - H*TAU)
 17				CONTINUE
C
C				CASE OF ROTATIONS Q < J <= N
C
				DO 18 J = IQ+1, N
				     G = A(IP,J)
				     H = A(IQ,J)
				     A(IP,J) = G - S*(H + G*TAU)
				     A(IQ,J) = H + S*(G - H*TAU)
 18				CONTINUE
				DO 19 J = 1,N
				     G = V(J,IP)
				     H = V(J,IQ)
				     V(J,IP) = G - S*(H + G*TAU)
				     V(J,IQ) = H + S*(G - H*TAU)
 19				CONTINUE
			    NROT = NROT + 1
                   END IF
 21		   CONTINUE
 22	    CONTINUE
C
C		UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
C
	  DO 23 IP = 1, N
	       B(IP) = B(IP) + Z(IP)
	       D(IP) = B(IP)
	       Z(IP) = 0.D0
 23	  CONTINUE
 24	CONTINUE
C
C		IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C		THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C		AND STOP
	if (imeg.eq.1) WRITE(UR7,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD
     1             NEVER HAPPEN'
	RETURN
	END
C **********************************************************************
		SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES D AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.
C
C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C ----------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
		K = I
		P = D(I)
		DO 11 J = I+1,N
		     IF(D(J) .GE. P) THEN
		          K = J
		          P = D(J)
		     END IF
 11		CONTINUE
		IF(K .NE. I) THEN
		   D(K) = D(I)
		   D(I) = P
		   DO 12 J = 1,N
		        P = V(J,I)
		        V(J,I) = V(J,K)
		        V(J,K) = P
 12	         CONTINUE
		END IF
 13	CONTINUE
	RETURN
	END
C ****************************************************************
C
C	THE FOLLOWING SUBROUTINES ARE SIMPLE UTILITY ROUTINES
C
C ********************************************************************
	SUBROUTINE ZEROV(V,SIZE)
C
C	THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C ---------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
	INTEGER SIZE
	REAL*8 V(0:SIZE-1)

	I = 0
10	IF (I .LE. SIZE-1) THEN
	   V(I) = 0.0
	   I = I+1
	GO TO 10
	END IF
	
	RETURN
	END
C **********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.0
C -------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 A(3,3)

        DO 1 I=1,3
        DO 1 J=1,3
             A(I,J) = 0.0
1	CONTINUE
	
	RETURN
	END
C **********************************************************************
      	SUBROUTINE ZEROMAT(A,N)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A N BY N MATRIX TO 0.0
C -------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 A(N,N)

        DO 1 I=1,N
        DO 1 J=1,N
             A(I,J) = 0.0
1	CONTINUE
	
	RETURN
	END
C ******************************************************************
	SUBROUTINE ONEM(A)
C
C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C ------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 A(3,3)

	DO 1 I = 1,3
	DO 1 J = 1,3
	     IF (I .EQ. J) THEN      
	         A(I,J)      = 1.0D0
	     ELSE
	         A(I,J)      = 0.0D0
	     END IF
1	CONTINUE

	RETURN
	END
C ******************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
C --
C --	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C --	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8  A(3,3), ATRANS(3,3)
C ------------------------------------------------------------------
C	COMPUTATION
C
        CALL ONEM(ATRANS)
	DO 1 I = 1, 3
	DO 1 J = 1, 3
		ATRANS(I,J) = A(J,I)
1	CONTINUE
	RETURN
	END
C ****************************************************************
	SUBROUTINE MPROD(A,B,C)
C --
C --	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C --	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3), B(3,3), C(3,3)
C ------------------------------------------------------------------
C	COMPUTATION
C
	DO 2 I = 1, 3
	DO 2 J = 1, 3
		C(I,J) = 0.0
		DO 1 K = 1, 3
			C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1		CONTINUE
2	CONTINUE
	RETURN
	END
C **************************************************************
	SUBROUTINE MPROD4(A,B,C)
C --
C --	THIS SUBROUTINE MULTIPLIES TWO 4 BY 4 MATRICES [A] AND [B],
C --	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(4,4), B(4,4), C(4,4)
C ------------------------------------------------------------------
C	COMPUTATION
C
	DO 2 I = 1, 4
	DO 2 J = 1, 4
		C(I,J) = 0.0
		DO 1 K = 1, 4
			C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1		CONTINUE
2	CONTINUE
	RETURN
	END
C ******************************************************************
	SUBROUTINE  DOTPM(A,B,C)
C
C	THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	SCALAR C.
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3),B(3,3),C

	C= 0.0
	DO 1 I = 1,3
	DO 1 J = 1,3
	     C = C + A(I,J)*B(I,J)
1	CONTINUE

	RETURN
	END
C********************************************************************
	SUBROUTINE  DOTSMV(A,B,C)

C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(6),B(6),C

	C= 0.0
	DO 1 I = 1,6
	     C = C + A(I)*B(I)
1	CONTINUE

C	DO 2 I = 4,6
C	     C = C + 2.*A(I)*B(I)
C2	CONTINUE

	RETURN
	END
C ******************************************************************
	SUBROUTINE MDET(A,DET)
C --
C --	THIS SUBROUTINE CALCULATES THE DETERMINANT
C --	OF A 3 BY 3 MATRIX [A].
C ------------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(3,3)
C ------------------------------------------------------------------
C	COMPUTATION
C
	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	+ A(1,2)*A(2,3)*A(3,1)
     +	+ A(1,3)*A(2,1)*A(3,2)
     +	- A(3,1)*A(2,2)*A(1,3)
     +	- A(3,2)*A(2,3)*A(1,1)
     +	- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END
C *******************************************************************
	SUBROUTINE M3INV(A,AINV)

C --	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX [A]
C --	AND PLACES THE RESULT IN [AINV]. 
C --	IF DET(A) IS ZERO, THE CALCULATION
C --	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C ----------------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)	
	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)
C
C		A(3,3)	-- THE MATRIX WHOSE INVERSE IS DESIRED.
C		DET		-- THE COMPUTED DETERMINANT OF [A].
C		ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C				   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C				   IS CALLED THE COFACTOR OF A(I,J).
C		AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C				   OBTAINED BY REPLACING EACH ELEMENT OF
C				   [A] BY ITS COFACTOR, AND THEN TAKING
C				   TRANSPOSE OF THE RESULTING MATRIX.
C		AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C				   [AINV] = [AADJ]/DET.
C ----------------------------------------------------------------------

	CALL MDET(A,DET)

	IF ( DET .EQ. 0.0 ) THEN
	    if (imeg.eq.1) WRITE(UR7,10)
	    CALL EXIT
	END IF

	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE

C ----------------------------------------------------------------------
C	FORMAT
C

 10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +       10X,'PROGRAM TERMINATED')
C----------------------------------------------------------------------
	RETURN
	END
C **********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
C --
C --	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C --	AND PLACES THE RESULT IN ACOFAC. 
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8  A(3,3), ACOFAC(3,3)
C ------------------------------------------------------------------
C	COMPUTATION
C
	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

C ****************************************************************
        SUBROUTINE INVAR(A,IA,IIA,IIIA)
C
C	THIS  SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	IA, IIA, IIIA OF	A  TENSOR [A].
C --------------------------------------------------------------------
C      	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA
C -------------------------------------------------------------------
C    	COMPUTATION
C
       DO 1 I=1,3
       DO 1 J=1,3
            AD(I,J) = A(I,J)
 1     CONTINUE

       IA = AD(1,1) + AD(2,2) + AD(3,3)
C
C			CALCULATE THE SQUARE OF [AD]
C
       CALL MPROD(AD,AD,AD2)
C
       IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )
C
       CALL  MDET(AD,DETA)

       IIIA = DETA
    
C
       RETURN
       END

C *********************************************************************
	SUBROUTINE TRACEM(A,TRA)
C
C	THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	AND STORES THE RESULT IN THE SCALAR TRA
C -------------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3),TRA

	TRA = A(1,1) + A(2,2) + A(3,3)

	RETURN 
	END
C ********************************************************************
	SUBROUTINE DEVM(A,ADEV)
C
C	THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	3 BY 3 MATRIX [A]
C ---------------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

	CALL TRACEM(A,TRA)
	CALL ONEM(IDEN)
	CALL ZEROM(ADEV)

	DO 1 I = 1,3
	DO 1 J = 1,3
	     ADEV(I,J) = A(I,J) - (1.0/3.0)*TRA*IDEN(I,J)
 1	CONTINUE

	RETURN
	END
C **********************************************************************
	SUBROUTINE EQUIVS(S,SB)
C
C	THIS SUBROUTINE CALCULATES THE EQUIVALENT TENSILE STRESS SB
C	CORESSPONDING TO A 3 BY 3 STRESS MATRIX [S]
C ---------------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

	SB = 0.0
	SDOTS = 0.0

	CALL DEVM(S,SDEV)
	CALL DOTPM(SDEV,SDEV,SDOTS)
	SB = DSQRT(1.5* SDOTS)

	RETURN
	END
C**************************************************************
	  SUBROUTINE PRTMAT(A,M,N)
	  IMPLICIT REAL*8(A-H,O-Z)
	  INTEGER M,N
	  REAL*8 A(M,N)	
        if (imeg.eq.1) then
	  DO 10 K=1,M
	    WRITE(UR7,'(2X,6F12.4,2X)') (A(K,L), L=1,N)
10      CONTINUE
        endif
        RETURN
        END


	  SUBROUTINE PRTJAC(A,M,N)
	  IMPLICIT REAL*8(A-H,O-Z)
	  INTEGER M,N
	  REAL*8 A(M,N)	
        if (imeg.eq.1) then
	  DO 10 K=1,M
	    WRITE(UR7,'(2X,6F12.3,2X)') (A(K,L), L=1,N)
10      CONTINUE
        endif
        RETURN
        END
C **********************************************************************
C
C	THE FOLLOWING SUBROUTINES ARE NEEDED FOR INTERFACING WITH 
C	ABAQUS
C
C **********************************************************************
          SUBROUTINE PUSHSV(SYM,VECT,IFLAG,NDI,NTENS)


C        IFLAG=1   CONVERTS A SYMMETRIC MATRIX WRITTEN AS A
C                  VECTOR VECT(6) TO A MATRIX SYM(3,3)
C        IFLAG=2   CONVERTS A SYMMETRIC MATRIX SYM(3,3) TO 
C                  THE CORRESPONDING VECTOR VECT(6)
C ----------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
         DIMENSION SYM(3,3),VECT(NTENS)
         NSHR=NTENS-NDI
         IF (IFLAG.EQ.1) THEN 
               CALL ZEROM(SYM)
               DO 15 I=1,NDI
15             SYM(I,I)=VECT(I)
               IF (NSHR.NE.0)    THEN
                   SYM(1,2)=VECT(NDI+1)
                   SYM(2,1)=VECT(NDI+1)
                        IF  (NSHR.NE.1)  THEN
                        SYM(1,3)=VECT(NDI+2)
                        SYM(3,1)=VECT(NDI+2)
                               IF  (NSHR.NE.2)  THEN
                               SYM(2,3)=VECT(NDI+3)
                               SYM(3,2)=VECT(NDI+3)
                               ENDIF
                        ENDIF
               ENDIF
        ELSE 
             IF (IFLAG.EQ.2) THEN 
                 DO 24 I=1,NTENS
24               VECT(I)=0.0
                 DO 25 I=1,NDI
25               VECT(I)=SYM(I,I)
                 IF (NSHR.NE.0) THEN
                       VECT(NDI+1)=SYM(1,2)
                       IF (NSHR.NE.1) THEN
                             VECT(NDI+2)=SYM(1,3)
                             IF (NSHR.NE.2) VECT(NDI+3)=SYM(2,3)
                       ENDIF
                 ENDIF
             ELSE
             WRITE(6,*) '** ERROR IN PUSHSV WRONG IFLAG  '
             ENDIF
         ENDIF

         RETURN
         END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	SUBROUTINE VECSMAT(V,A)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION A(3,3),V(6)
	DO 10 I = 1,3
10	  A(I,I) = V(I)
	A(1,2) = V(4)  
        A(2,1) = V(4)
	A(1,3) = V(5) 
        A(3,1) = V(5)
	A(2,3) = V(6) 
        A(3,2) = V(6)
	RETURN
	END
C 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C 
      SUBROUTINE SMATVEC(A,V)
 	IMPLICIT REAL*8(A-H,O-Z)       
      DIMENSION A(3,3),V(6)
      DO 10 I = 1,3
10      V(I) = A(I,I)
        V(4) = A(1,2)
        V(5) = A(1,3)
        V(6) = A(2,3)
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     	SUBROUTINE ZEROMOD(A)

	IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(3,3,3,3)

      DO 1 I=1,3
      DO 1 J=1,3
	DO 1 K=1,3
	DO 1 L=1,3
             A(I,J,K,L) = 0.0
 1	CONTINUE
	
	RETURN
	END	
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	SUBROUTINE EULANG(Q,N,TH,PHI,OM,ICRYS,IEULERERR)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION Q(N,N)
c-twq        PI = 4.0*DATAN(1.0D0)
	ICHK=0
        IEULERERR = 0
	IF(DABS(Q(3,3))-1.0.GT.1.0D-6)THEN
          if (imeg.eq.1)then
           WRITE(UR7,*)'Q',((Q(I,J),J=1,3),I=1,3)
	   WRITE(UR7,*) 'Q(3,3) > 1.0'
         endif
	   IEULERERR = 1
	    RETURN
	ENDIF
	DO 10 I = 1,3
	DO 10 J = 1,3
10	  IF(DABS(Q(I,J)).LT.1.0D-6)Q(I,J)=0.0
	IF(DABS(DABS(Q(3,3))-1.0).LT.1.0D-6)THEN
          CALL EULCHECK1(Q,3,TH,PHI,OM,ICHK)
	  IF(ICHK.NE.1)GO TO 20
	  RETURN
	ENDIF
	TH = DACOS(Q(3,3))
	STH = DSIN(TH)
	OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
	PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
	CALL EULCHECK(Q,3,TH,PHI,OM,ICHK)
	IF(ICHK.EQ.1)RETURN
	TH = 2.0*PI-TH
	STH = DSIN(TH)
	OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
	PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
	CALL EULCHECK(Q,3,TH,PHI,OM,ICHK)
20	IF(ICHK.NE.1)THEN
          if (imeg.eq.1)then
	    WRITE(9,*)'ICRYS = ',ICRYS
	    WRITE(9,*)'Q',((Q(J,K),K=1,3),J=1,3)
            WRITE(UR7,*) 'FAILED TO FIND EULER ANGLES'
          endif
	    IEULERERR = 1
	    RETURN
	ENDIF
	RETURN
	END
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE EULCHECK(Q,N,TH,PHI,OM,ICHK)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION Q(N,N)
	TOL=1.0D-3
	A = DCOS(PHI)*DCOS(OM)-DSIN(PHI)*DSIN(OM)*DCOS(TH)
	B = -DSIN(OM)*DCOS(PHI)-DCOS(OM)*DSIN(PHI)*DCOS(TH)
	C = DCOS(OM)*DSIN(PHI)+DSIN(OM)*DCOS(PHI)*DCOS(TH)
	D = -DSIN(PHI)*DSIN(OM)+DCOS(PHI)*DCOS(OM)*DCOS(TH)
        IF(DABS(A-Q(1,1)).LT.TOL.AND.DABS(B-Q(2,1)).LT.TOL.
     +    AND.DABS(C-Q(1,2)).LT.TOL.AND.DABS(D-Q(2,2)).LT.TOL)ICHK=1
	RETURN
	END
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE EULCHECK1(Q,N,TH,PHI,OM,ICHK)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION Q(N,N)
	TOL=1.0D-3
	Q(3,3) = 1.0*Q(3,3)/DABS(Q(3,3))
	TH = DACOS(Q(3,3))
	IF(DABS(Q(1,3)).GT.TOL)RETURN
	IF(DABS(Q(2,3)).GT.TOL)RETURN
	IF(DABS(Q(3,1)).GT.TOL)RETURN
	IF(DABS(Q(3,2)).GT.TOL)RETURN
	IF(Q(3,3).EQ.1.0)THEN
	   IF(DABS(Q(1,1)-Q(2,2)).GT.TOL)RETURN
	   IF(DABS(Q(1,2)+Q(2,1)).GT.TOL)RETURN
	ELSE
	   IF(DABS(Q(1,1)+Q(2,2)).GT.TOL)RETURN
	   IF(DABS(Q(1,2)-Q(2,1)).GT.TOL)RETURN
	ENDIF
	PHI = DATAN2(Q(1,2),Q(1,1))
	OM = 0.0
	ICHK = 1
	RETURN
	END
c
C ****************************************************************
	SUBROUTINE VPROD(A,B,C)
C 
C 	THIS SUBROUTINE MULTIPLIES 3 BY 3 MATRICES [A] AND 3 BY 1 VECTOR [B],
C 	AND PLACE THEIR PRODUCT IN VECT [C]. 
C -----------------------------------------------------------------
C	VARIABLES
C
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3), B(3), C(3)

	DO I=1,3
		C(I)=0.
		DO J=1,3
			C(I)=C(I)+A(I,J)*B(J)
		ENDDO
	ENDDO

	RETURN
	END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	
	SUBROUTINE DEG2RAD(TH,ITAG)
c
	IMPLICIT REAL*8(A-H,O-Z)
C
C	DEGREE <-> RADIAN
	
      PI0=2.*ASIN(1.0)
	
	IF(ITAG.EQ.1) THEN

		TH=TH*PI0/180.0

	ELSE 

		TH=TH*180./PI0
	ENDIF

	RETURN
	END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	subroutine uptex(q,theta,phi,omega,itex_err)
c
	implicit real*8(a-h,o-z)
	dimension q(3,3)
     
c-twq	pi = 4.0*datan(1.0d0)

c	First check Q(3,3): cos(theta)	
	


c	CASE A : write error

	if((Q(3,3)-1.d0).gt.1.d-6) then

		itex_err=1

          if (itex.eq.1)	write(UR3,*) 'Check error in uptex, Q33 is 
     1                                      greater than 1.'
		return

	elseif((Q(3,3)+1.d0).lt.-1.d-6) then ! error

		itex_err=1
	if (itex.eq.1) write(UR3,*) 'Check error in uptex, Q33 is less than -1.'
		return

	
c	CASE B : cos(theta)=1 or -1	
	
	elseif(dabs(Q(3,3)-1.d0).lt.1.d-6.or.
	1		dabs(Q(3,3)+1.d0).lt.1.d-6) then 

		do i=1,3
		do j=1,3
			if(dabs(Q(i,j)).lt.1.d-6) Q(i,j)=0.d0
		enddo
		enddo

		if(Q(3,3).gt.0.d0) then
			Q(3,3)=1.d0
		elseif(Q(3,3).lt.0.d0) then
			Q(3,3)=-1.d0
		endif

		theta=dacos(Q(3,3)) ! 0. or pi
		
		!Double check, check error

		if(dabs(Q(1,3)).gt.1.d-3.or.dabs(Q(2,3)).gt.1.d-3.or.
	1		dabs(Q(3,1)).gt.1.d-3.or.dabs(Q(3,2)).gt.1.d-3) then

		itex_err=1
		if (itex.eq.1) write(UR3,*) 'Check error in uptex, when cos(theta)=1
     1                                or -1'
		return

		endif

		if(Q(3,3).eq.1.d0) then

			if(dabs(Q(1,1)-Q(2,2)).gt.1.d-3.or.
	1			dabs(Q(2,1)+Q(1,2)).gt.1.d-3) then

				itex_err=1
				if (itex.eq.1) write(UR3,*) 'Check error in uptex,when cos(th)=1'
				return

			endif
		
		elseif(Q(3,3).eq.-1.d0) then

			if(dabs(Q(1,1)+Q(2,2)).gt.1.d-3.or.
	1		   dabs(Q(2,1)-Q(1,2)).gt.1.d-3) then

				itex_err=1
				if (itex.eq.1) write(UR3,*) 'Check error in uptex,when cos(th)=-1'
				return

			endif

		endif

		omega=0.d0
		phi=datan2(Q(1,2),Q(1,1))

		itex_err=0 ! pass
		return	



c	CASE C : General case


	else ! general case
		
c		write(UR3,*) 'q33',Q(3,3)
		
		theta=dacos(Q(3,3))
		st=dsin(theta)
		omega=datan2((Q(1,3)/st),(Q(2,3)/st))
		phi=  datan2((Q(3,1)/st),(-Q(3,2)/st))

c		write(UR3,*)'check',theta,omega,phi
	
		! Double check
		
		Q11=dcos(phi)*dcos(omega)-dsin(phi)*dsin(omega)*dcos(theta)
		Q12=dsin(phi)*dcos(omega)+dcos(phi)*dsin(omega)*dcos(theta)
		Q21=-dcos(phi)*dsin(omega)-dsin(phi)*dcos(omega)*dcos(theta)
		Q22=-dsin(phi)*dsin(omega)+dcos(phi)*dcos(omega)*dcos(theta)

		if(dabs(Q(1,1)-Q11).gt.1.d-3.or.dabs(Q(2,2)-Q22).gt.1.d-3.or.
	1	dabs(Q(1,2)-Q12).gt.1.d-3.or.dabs(Q(2,1)-Q21).gt.1.d-3) then

			theta=2.d0*pi-theta
			st=dsin(theta)
			omega=datan2((Q(1,3)/st),(Q(2,3)/st))
			phi=  datan2((Q(3,1)/st),(-Q(3,2)/st))

		endif

		Q11=dcos(phi)*dcos(omega)-dsin(phi)*dsin(omega)*dcos(theta)
		Q12=dsin(phi)*dcos(omega)+dcos(phi)*dsin(omega)*dcos(theta)
		Q21=-dcos(phi)*dsin(omega)-dsin(phi)*dcos(omega)*dcos(theta)
		Q22=-dsin(phi)*dsin(omega)+dcos(phi)*dcos(omega)*dcos(theta)

		if(dabs(Q(1,1)-Q11).gt.1.d-3.or.dabs(Q(2,2)-Q22).gt.1.d-3.or.
	1	dabs(Q(1,2)-Q12).gt.1.d-3.or.dabs(Q(2,1)-Q21).gt.1.d-3) then

			itex_err=1
			if (itex.eq.1) write(UR3,*) 'Check error in uptex, general case'
			return

		endif

		itex_err=0 !pass

		return


	endif

	return
	end

c======================================================================
	subroutine inout(ce,nce,ntgt,crd,nnte,itag)
      use modname	
	IMPLICIT REAL*8(A-H,O-Z)

!	Include 'D:\temp-wagoner102\COMMON'

	dimension ce(MAXNOEL,8,3)
	dimension nce(ntgt)
	dimension amm(3,3),av(3,3),amt(3,3)
	dimension a(3),b(3),c(3),crd(3)
	dimension x(3),s(3)
		
	itag=1

	do k=1,ntgt
		nelem=nce(k)

		do i=1,3
			amm(1,i)=ce(nelem,2,i)-ce(nelem,1,i)
			amm(2,i)=ce(nelem,4,i)-ce(nelem,1,i)
			amm(3,i)=ce(nelem,5,i)-ce(nelem,1,i)
		enddo
		
		do i=1,3	
			x(i)=crd(i)-ce(nelem,1,i)
		enddo

		
		do i=1,3
			do j=1,3
				amt(i,j)=amm(j,i) !!
			enddo
		enddo

		call m3inv(amt,av)
		call VPROD(av,x,s)
	
		if(s(1).lt.0..or.s(2).lt.0..or.s(3).lt.0.) then
			itag=1
			goto 100
		endif
c		
		do i=1,3
			amm(1,i)=ce(nelem,8,i)-ce(nelem,7,i)
			amm(2,i)=ce(nelem,6,i)-ce(nelem,7,i)
			amm(3,i)=ce(nelem,3,i)-ce(nelem,7,i)
		enddo
		
		do i=1,3	
			x(i)=crd(i)-ce(nelem,7,i)
		enddo

		
		
		do i=1,3
			do j=1,3
				amt(i,j)=amm(j,i) !!
			enddo
		enddo
		
		call m3inv(amt,av)
		call VPROD(av,x,s)
		
		if(s(1).lt.0..or.s(2).lt.0..or.s(3).lt.0.) then
			itag=1
			goto 100
		else
			itag=0
			nnte=nelem   ! find target
			return
		endif

100	continue
	enddo

	return
	end subroutine
c
c
c-----------------------------------------------------------------------
	subroutine get_rss(mslip,nnte,rss)
      use modname
	IMPLICIT REAL*8(A-H,O-Z)

!	Include 'D:\temp-wagoner102\COMMON'
	
		rss=0.
		do ipt=1,maxnpt
			rss=rss+rshear(nnte,ipt,mslip)
		enddo
		rss=rss/maxnpt

	return		

	end subroutine

c
C---------------------------------------------------------------------------------
C--------FINDING NEIGHBORING ELEMENTS AND KEEP THE ELEMENT NUMBER-----------------
c---------------------------------------------------------------------------------
	subroutine elcandi(NSTORE1,nelc1)
      use modname
	implicit real*8(a-h,o-z)

!	Include 'D:\temp-wagoner102\COMMON'
c	
	parameter(nlayer=2)
	integer nel(maxnode,maxkount),nkount(maxnode)
	integer nelc1(maxnoel,8) 
c-twq      
c-twq      integer,allocatable:: nelc1(:,:)     
	integer NSTORE1(MAXNOEL,500)
C      DIMENSION NK1(MAXNOEL)	
c
c	STEP 1: MAKE A FILE WHICH CONTAINS THE ELEMENTS CONNECTED WITH NODES
c
c-twq	open(61, name='D:\temp-wagoner102\TESTA.ele')
      read(UR2,*)
c	read(UR0,*) nelem !NUMBER OF ELEMENTS
	nelem=MAXNOEL

	do iel=1,nelem
		read(UR2,*) ntemp, (nelc1(iel,j),j=1,8)
c          if (itex.eq.1)then
c		if(ntemp.ne.iel) write(*,*) 'Error 0'
c          endif
c     write 'error 0' if 'ntemp' is not equal to 'iel'
	enddo
	
	do inode=1,maxnode
		kount=0
		do jelem=1,nelem
			do k=1,8
				mnode=nelc1(jelem,k)
				if(mnode.eq.inode) then
					kount=kount+1
					nel(inode,kount)=jelem
				endif
			enddo
		enddo

		if(kount.eq.0) goto 100
		
		nkount(inode)=kount
c-twq			rewind(61)
c-twq		read(61,*) nelem
	enddo

100	continue

c	STEP 2: SEARCH NEIGHBORING ELEMENTS

	do kel=1,nelem
	
C	write(UR10,*)kel

		nstore1(kel,1)=kel ! lst stored element is always target element itself
		kount=1
		do ilayer=1,nlayer
			if(ilayer.eq.1) then
			
				do inode=1,8
					nn=nelc1(kel,inode) ! connected node
					nkk=nkount(nn)
					do kk=1,nkk ! # of connected element
						itag=0
						do kkk=1,kount
							if(nel(nn,kk).eq.nstore1(kel,kkk)) itag=1
						enddo
						if(itag.eq.0) then
							kount=kount+1
							nstore1(kel,kount)=nel(nn,kk)
						endif
					enddo
				enddo
				
				kount_e=kount
			
			else

				do is=2,kount_e
					nn=nstore1(kel,is)
					do inode=1,8
						mm=nelc1(nn,inode)
						nkk=nkount(mm)
						do kk=1,nkk
							itag=0
							do kkk=1,kount
								if(nel(mm,kk).eq.nstore1(kel,kkk)) 
	1							itag=1
							enddo
							if(itag.eq.0) then
								kount=kount+1
								nstore1(kel,kount)=nel(mm,kk)
							endif
						enddo
					enddo
				enddo
			endif
				
		enddo

		nk(kel)=kount
				
	enddo
	
	return
	end subroutine

c---------------------------------------------------------------------------------
c-------READS ELEMENTS FOR EACH GRAIN---------------------------------------------
C---------------------------------------------------------------------------------
C     THIS SUBROUTINE SAVES START AND END ELEMENTS (IN ABAQUS) FOR EACH GRAIN 
C

	
c#####################################################################
	SUBROUTINE REDIST3(NGD,DENS0g,BF,leng1,leng2,leng3,GAML,rsst,xc,
     1 islip)	
            use modname
!	Include 'E:\SDM\SDM\common'
		
	INTEGER igrid,ngrid,ngd
	INTEGER i,j,k,ii,l,islip

	REAL(8) bv,xnu,g,xk,sum,dxq,thk,rsst,tauobs
	REAL(8) dens1(maxkount),bf(maxkount),rrss(maxkount)
	real(8) dens0g(maxkount),dent(maxkount),dentt(maxkount)

	REAL(8) leng1(maxkount),leng2(maxkount),leng3(maxkount)
	REAL(8) f(maxkount),sig(maxkount),rrt(maxkount),GAML(maxkount)
	REAL(8) r1(maxkount,maxkount)
	REAL(8) gtot(maxkount,maxkount),g11,g12
      
      ! added by jisheng
      real(8) xc(MAXNOEL,3),vectorg(3),vectorl(3),
     1  D0(MAXSLIP,3),rot1(3,3),r1c,r2c,r3c

	

	xk=2.*pi*(1.-0.3d0)/bs/shm
	
	rsst=abs(rsst)
	ii=0

C----------------REDISTRIBUTION OF DISLOCATIONS2---------
	
	dent(1:ngd)=0 !number of mobile dislocations 
	dentt(1:ngd)=0
c	
c	!Calculate number of mobile dislocation 
	
	IF (NGD.EQ.1) THEN
		DENTT(1)=LENG2(I)/bs*GAML(I)
	ELSE

		DO I=1,NGD-1
		dent(i)=(leng2(i)+leng2(i+1))/2/bs
	1	*(gaml(i)+gaml(i+1))/2
		ENDDO

		dentt(1)=dent(1)
		dentt(ngd)=0	!last elements
c	
		do i=2,ngd-1
		dentt(i)=dent(i)-dent(i-1)
		enddo

	ENDIF

	DO i=1,ngd	
		dens0g(i)=dens0g(i)+dentt(i)
	ENDDO

c----------------BACK STRESS CALCULATION: FINITE LINE LENGTH-------------
	
	f(1:(ngd))=0.
	rrt(1:ngd)=0

	IF (ngd.eq.1) THEN
		bf(1)=0
	ELSE
		
		rrt(1)=0.5*leng1(1)
		DO i=2,ngd
		rrt(i)=rrt(i-1)+0.5*leng1(i)+0.5*leng1(i-1)
		ENDDO
	
		do i=1,ngd

		sum=0
			do j=1,ngd
				if (i.ne.j) then	
      !!!!!!!!jisheng!!!!!!!
                  vectorg(1:3)=xc(nte(j),1:3)-xc(nte(i),1:3)



          do k = 1,3   
          do l = 1,3 		             
			 rot1(k,l)= QMAT(j,NPT1,k,l)  
          enddo
          enddo   

           
		do k=1,3
			vectorl(k)=0.
		do l=1,3
                     
		vectorl(k)=vectorl(k)+rot1(l,k)*vectorg(l) 
          ! rot1 is transpose of QMAT, QMAT rotates coordinates from local to global
		enddo
          enddo

      
      r1c=B0(islip,1)*vectorl(1)+B0(islip,2)*vectorl(2) !the distance calculated from centrol coordinates
     1 +B0(islip,3)*vectorl(3)
      r2c=C0(islip,1)*vectorl(1)+C0(islip,2)*vectorl(2)
     1 +C0(islip,3)*vectorl(3) 
      call crospd(B0(islip,1:3),C0(islip,1:3),D0(islip,1:3))
      r3c=D0(islip,1)*vectorl(1)+D0(islip,2)*vectorl(2)
     1 +D0(islip,3)*vectorl(3) 

      !!!!jisheng!!!!!!
                        
				r1(i,j)=(rrt(i)-rrt(j))
				
				g11=(0.5*leng3(i)-0.5*leng3(j))**2
				g12=(0.5*leng3(i)+0.5*leng3(j))**2
				gtot(i,j)=sqrt(g11+r1(i,j)**2)-sqrt(g12+r1(i,j)**2)				
				gtot(i,j)=-2*gtot(i,j)
!deleted jisheng  f(i)=dens0g(j)/xk/2/R1(I,J)*gtot(i,j)/leng3(i)
                     
                     f(i)=dens0g(j)/xk/2/r1c*gtot(i,j)/leng3(i)

!                     f(i)=-sign(f(i),dens0g(i))

                      else
				f(i)=0
				endif	
			sum=sum+f(i)	
			enddo
		bf(i)=sum

		enddo

	endif

C-----------------------------------------------------------------------
803	continue


	return
	end subroutine	
c****************************************************************      
c-twq
      SUBROUTINE transmission(htrans1,elem,egrain,gangle2,g,plane)
         use modname 
c
      integer itemp,i,i1,i2,j,k
      integer nid(4),eid1,eid2
c-twq      real(8),allocatable:: plane(:,:),g(:,:)
      real(8) p1,p,p2,gbn(2),q1,q,q2,ang1,ang,ang2
      real(8),allocatable:: node(:,:),gbnormal(:,:,:) 
      integer,allocatable:: neigh_elem(:,:),
     1    nneigh(:)
      real(8) n,nmax,g1(3),gi(3),l1(3),li(3),glocal(3)
      integer elemmax,slipmax,elem1,slip1,islip,ielem,eid
      real(8) plocal(3),pglobal1(3),pglobali(3)
      integer iline,ibcount
      real(8) lmax,gmax
      REAL(8) G1MAX(3),GIMAX(3),P1MAX(3),PIMAX(3)
c      
      real(8) htrans1(maxnoel,MAXSLIP)
      integer egrain(maxnoel)
      real(8) gangle2(MAXGRAIN,3)
      integer elem(MAXNOEL,8)
      real(8) plane(MAXSLIP,3),g(MAXSLIP,3)
      
      allocate(node(maxnode,3))      
      allocate(gbnormal(MAXNOEL,4,2))
      allocate(nneigh(MAXNOEL)); nneigh=0
      allocate(neigh_elem(MAXNOEL,4))    
      
c-twq      integer,allocatable::nb1(:,:),nbind1(:)      
c-twq      allocate(node(maxnode,3))
c-twq      allocate(elem(MAXNOEL,8))
c-twq      allocate(gangle2(MAXGRAIN,3))
c-twq      allocate(egrain(MAXNOEL))
c-twq      allocate(gbnormal(MAXNOEL,4,2))
c-twq      allocate(nneigh(MAXNOEL)); nneigh=0
c-twq      allocate(neigh_elem(MAXNOEL,4))
c-twq      allocate(plane(MAXSLIP,3))
c-twq      allocate(g(MAXSLIP,3))
c-twq      allocate(nb(MAXGRAIN,MAXNOEL))
c-twq      allocate(nbind(MAXGRAIN))

      ! read node
c-twq      open(87,name='data.dat')
c-twq      open(UR10,name='testa1.nod')

      read(UR2,*)

      do i=1,maxnode      
         read(UR2,*) itemp,node(itemp,1),node(itemp,2),node(itemp,3)
      enddo
      
      open(UR11,file=Trim(AdjustL(FILEroutine))//'transmissivity.out')
      
      ! build gb neighbor list
      do eid1=1,MAXNOEL
      do eid2=eid1+1,MAXNOEL
      ncount=0
      do i=1,8
      do j=1,8
      if(egrain(eid1).ne.egrain(eid2)) then
      if(elem(eid1,i).eq.elem(eid2,j)) then
      ncount=ncount+1
      nid(ncount)=elem(eid1,i)
      if(ncount.eq.4)then
      goto 9515
      endif
      endif

      endif
      enddo
      enddo
      enddo !eid2

9515  continue

      if(ncount.eq.4) then ! neighboring
      gbn(1)=node(nid(1),2)-node(nid(2),2)
      gbn(2)=node(nid(2),1)-node(nid(1),1)

      if(abs(gbn(1)).le.1.0e-10)then
      gbn(1)=node(nid(1),2)-node(nid(3),2)
      gbn(2)=node(nid(3),1)-node(nid(1),1)			
      endif

      if(abs(gbn(1)).le.1.0e-10)then
      gbn(1)=node(nid(1),2)-node(nid(4),2)
      gbn(2)=node(nid(4),1)-node(nid(1),1)			
      endif			

      nneigh(eid1)=nneigh(eid1)+1
      neigh_elem(eid1,nneigh(eid1))=eid2
      gbnormal(eid1,nneigh(eid1),1)=gbn(1)
      gbnormal(eid1,nneigh(eid1),2)=gbn(2)

      nneigh(eid2)=nneigh(eid2)+1
      neigh_elem(eid2,nneigh(eid2))=eid1
      gbnormal(eid2,nneigh(eid2),1)=gbn(1)
      gbnormal(eid2,nneigh(eid2),2)=gbn(2)
      write(*,*) 'neighbor',eid1,eid2
      endif
      enddo !eid1

      !! input
      elem1=1
      slip1=1

      do elem1=1,MAXNOEL
      ang1=gangle2(egrain(elem1),1)/180.d0*3.141592
      ang=gangle2(egrain(elem1),2)/180.d0*3.141592
      ang2=gangle2(egrain(elem1),3)/180.d0*3.141592
      do slip1=1,MAXSLIP

      nmax=0
      !elemmax
      !slipmax

      glocal=g(slip1,:)
      call local2global(ang1,ang,ang2,glocal,g1)
      call normalize(g1)

      plocal=plane(slip1,:)

      call local2global(ang1,ang,ang2,plocal,pglobal1)
      call normalize(pglobal1)

      do ineigh=1,nneigh(elem1)
      if(pglobal1(3).ne.0) then
      l1(1)=gbnormal(elem1,ineigh,2)
      l1(2)=-gbnormal(elem1,ineigh,1)
      l1(3)=-(pglobal1(2)*l1(2)+pglobal1(1)*l1(1))/pglobal1(3)
      call normalize(l1)
      else
      l1(1)=0;l1(2)=0;l1(3)=1.d0
      endif

      eid=neigh_elem(elem1,ineigh)
      q1=gangle2(egrain(eid),1)/180.d0*3.141592
      q=gangle2(egrain(eid),2)/180.d0*3.141592
      q2=gangle2(egrain(eid),3)/180.d0*3.141592
      do islip=1,MAXSLIP
      glocal=g(islip,:)

      call local2global(q1,q,q2,glocal,gi)
      call normalize(gi)
      plocal=plane(islip,:)

      call local2global(q1,q,q2,plocal,pglobali)
      call normalize(pglobali)

      if(pglobali(3).ne.0) then
      li(1)=l1(1)
      li(2)=l1(2)
      li(3)=-(pglobali(2)*li(2)+pglobali(1)*li(1))/pglobali(3)
      call normalize(li)
      else
      li(1)=0;li(2)=0;li(3)=1.d0
      endif

      n=dabs( (l1(1)*li(1)+l1(2)*li(2)+l1(3)*li(3))*(g1(1)*gi(1)
     1    +g1(2)*gi(2)+g1(3)*gi(3)) )

      if(n.gt.nmax) then
      nmax=n
      elemmax=eid
      slipmax=islip
      lmax=lmulti
      gmax=gmulti
      G1MAX(1)=G1(1)
      G1MAX(2)=G1(2)
      G1MAX(3)=G1(3)
      GIMAX(1)=GI(1)
      GIMAX(2)=GI(2)
      GIMAX(3)=GI(3)	
      P1MAX(1)=PGLOBAL1(1)
      P1MAX(2)=PGLOBAL1(2)
      P1MAX(3)=PGLOBAL1(3)
      PIMAX(1)=PGLOBALI(1)
      PIMAX(2)=PGLOBALI(2)
      PIMAX(3)=PGLOBALI(3)			
      endif
      enddo
      enddo


      IF (NMAX.NE.0) THEN
c-twq      WRITE(16,*) ELEM1
 
      write(UR11,*) elem1,slip1,nmax
c-twq
      htrans1(elem1,slip1)=nmax
      ELSE
      write(UR11,*) elem1,slip1,1
c-twq
      htrans1(elem1,slip1)=1
      ENDIF
      
 

      enddo

      enddo

c      close(87)

	CLOSE(UR11)
	return
	end subroutine	


      !----!----!----!----!----!----!----!----!----!----!----
      subroutine local2global(pp1,pp,pp2,vlocal,vglobal)
      real(8) pp1,pp,pp2
      real(8) vlocal(3),vglobal(3)
      real(8) p1mat(3,3)
      !c-twq 
c      parameter (pi=3.141592654)

      !p1mat=0; pmat=0; p2mat=0
      p1mat(1,1)=0
      p1mat(1,2)=0
      p1mat(1,3)=0
      p1mat(2,1)=0
      p1mat(2,2)=0
      p1mat(2,3)=0
      p1mat(3,1)=0
      p1mat(3,2)=0
      p1mat(3,3)=0
c
      p1mat(1,1)=DCOS(pp1)*DCOS(pp2)-DSIN(pp1)*DCOS(pp)*DSIN(pp2)
      p1mat(1,2)=-DCOS(pp1)*DSIN(pp2)-DSIN(pp1)*DCOS(pp)*DCOS(pp2)
      p1mat(1,3)=DSIN(pp1)*DSIN(pp)
      p1mat(2,1)=DSIN(pp1)*DCOS(pp2)+DCOS(pp1)*DCOS(pp)*DSIN(pp2)
      p1mat(2,2)=-DSIN(pp1)*DSIN(pp2)+DCOS(pp1)*DCOS(pp)*DCOS(pp2)
      p1mat(2,3)=-DCOS(pp1)*DSIN(pp)
      p1mat(3,1)=DSIN(pp)*DSIN(pp2)
      p1mat(3,2)=DSIN(pp)*DCOS(pp2)
      p1mat(3,3)=DCOS(pp)

      vglobal(1)=p1mat(1,1)*vlocal(1)+p1mat(1,2)*vlocal(2)
     1    +p1mat(1,3)*vlocal(3)
      vglobal(2)=p1mat(2,1)*vlocal(1)+p1mat(2,2)*vlocal(2)
     1    +p1mat(2,3)*vlocal(3)
      vglobal(3)=p1mat(3,1)*vlocal(1)+p1mat(3,2)*vlocal(2)
     1    +p1mat(3,3)*vlocal(3)
c
      return
      end subroutine
      !----!----!----!----!----!----!----!----!----!----!----
c
      subroutine normalize(v)
      real(8) v(3),vn
      vn=dsqrt(v(1)**2+v(2)**2+v(3)**2)
      v=v/vn
      return
      end subroutine
!************************************************************************
!twq   read form abaqus.inp  
      SUBROUTINE inp_info() 
      use modname
      IMPLICIT REAL*8(A-H,O-Z)    
      character(79) keywd,keywd0
      integer status,i,j,iel,itemp,ntemp
      integer ig,is,ie
      integer NyNode,NyEle,NyGB,NyGrain
C      integer maxnode,MAXNOEL
c      parameter(maxnode=450,MAXNOEL=196)     
      real(8) node(maxnode,3)
      integer nelc1(MAXNOEL,8)
c      
      NyNode=0
      NyEle=0
      NyGB=0
      NyGrain=0
c   
!     READS from anaqus.inp     
      READ(UR0,*)
	READ(UR0,'(A)') FILEINP    ! inp	
c
C      OPEN(UR1,file='D:\temp-wagoner102\fourgrains-sjtu.inp')
	OPEN(UR1,name=Trim(AdjustL(FILEroutine))//Trim(AdjustL(FILEINP)))
C	open(UR2,name='D:\temp-wagoner102\abaqus_info.dat')   !output of inp 
	open(UR2,file=Trim(AdjustL(FILEroutine))//'abaqus_info.dat')   !output of inp 
c 
! Grain inf.      
      do while (.true.)
      read(UR1,'(A79)',iostat=status) keywd
      if(status/=0) exit
      if(NyGrain==0)then
      if(keywd(1:19)=='*Elset, elset=GRAIN') then         
           NyGrain=1
           write(UR2,*)'*Grain++++++++++++++++++++++++++++++++++++++++++
     1+++++++++++++++++++++++++'           
           write(UR2,102)'grain #','st ele', 'end ele' 
		 backspace(UR1)
           do i=1,MAXGRAIN 
		   read(UR1,'(A79)') keywd0
		   read(UR1,*) is,ie
             write(UR2,103) keywd0(20:20),is,ie             
	     enddo        
         endif
      endif
	enddo

! Element inf. 
	rewind(UR1)
      do while (.true.)
	read(UR1,'(A79)',iostat=status) keywd
      if(status/=0) exit
      if(keywd(1:8)=='*Element') then
        if(NyEle==0)then
         NyEle=1
         write(UR2,*)'*Element++++++++++++++++++++++++++++++++++++++++++
     1+++++++++++++++++++++++++'         
	   do iel=1,MAXNOEL
		  read(UR1,*) ntemp, (nelc1(iel,j),j=1,8)
		  if(ntemp.ne.iel) write(*,*) 'Error element'
            write(UR2,101)ntemp,(nelc1(iel,j),j=1,8)
	   enddo
        endif
      endif
	enddo

! Node inf.     
	rewind(UR1)
      do while (.true.)
	read(UR1,'(A79)',iostat=status) keywd
      if(status/=0) exit
      if(NyNode==0)then
      if(keywd(1:5)=='*Node') then        
           NyNode=1
           write(UR2,*)'*node++++++++++++++++++++++++++++++++++++++++++
     1+++++++++++++++++++++++++'          
           do i=1,maxnode      
            read(UR1,*) itemp,node(itemp,1),node(itemp,2),node(itemp,3)
		  if(itemp.ne.i) write(*,*) 'Error node'
           write(UR2,100)itemp,node(itemp,1),node(itemp,2),node(itemp,3)
           enddo
        endif
      endif
	enddo

! GB inf.   
	rewind(UR1)
      do while (.true.)
	read(UR1,'(A79)',iostat=status) keywd
      if(status/=0) exit
      if(keywd(1:17)=='*Elset, elset=GB1') then
         if(NyGB==0)then
            NyGB=1
            write(UR2,*)'*GB++++++++++++++++++++++++++++++++++++++++++
     1+++++++++++++++++++++++++'            
            backspace(UR1)
            do 
		    read(UR1,"(A79)")keywd
              if(keywd(1:20)=='*Elset, elset=GRAIN1') goto 204
                 write(UR2,*)keywd
	      enddo        
         endif
      endif
      enddo   !enddo while(.true.)

! write format
100   format(i5,3f18.9) 
101   format(9i8) 
102   format(3a12)  
103   format(a12,2i12)       
c      

c      
204      close(UR1)
	   close(UR2)
c  
c      
	RETURN
	END
