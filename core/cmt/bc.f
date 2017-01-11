C> @file Routines for boundary conditions
      subroutine bcflux(flux,agradu,qminus)
! Need proper indexing and nekasgn & cmtasgn calls
      include 'SIZE'
      include 'INPUT'
      include 'DG'
!     include 'NEKUSE'
      include 'TSTEP' ! wait how do we know what ifield is?
      integer e,eq,f
      real flux  (nx1*nz1,2*ndim,nelt,toteq)
      real agradu(nx1*nz1,2*ndim,nelt,toteq)
!     real qminus(nx1*nz1,2*ndim,nelt,nqq) ! include CMTDATA?
      real qminus(*) ! 'scuse me. comin' through
      common /nekcb/ cb
      character*3 cb

      nfaces=2*ndim
      nxz=nx1*nz1
      ifield=1

      do e=1,nelt
         do f=1,nfaces
            if (cbc(f, e, ifield).ne.'E  '.and.
     >          cbc(f, e, ifield).ne.'P  ') then ! cbc bndy
               cb=cbc(f,e,ifield)
               if (cb .eq. 'I  ') then ! NEUMANN CONDITIONS GO HERE
!-------------------------------------------------------------
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
                  call rzero(flux(1,f,e,1),nxz)
                  do eq=2,ndim+1
                     call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
                  enddo
! METHOD "B", ADIABATIC NO-SLIP
                  call rzero(flux(1,f,e,toteq),nxz)
! METHOD "A", ADIABATIC NO-SLIP augments with viscous work. triage below
! because, predictably, NOW I need to computate AgradU on surfaces and I don't
! have general code for that.
                  call a5adiabatic_wall(flux(1,1,1,toteq),f,e,agradu,
     >                                  qminus)
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
!                 cbu=cb
!                 do eq=1,toteq
!                    call userflux(flux(1,f,e,eq)) ! replace this with userbc
!                 enddo
               elseif (cb .eq. 'SYM') then
                  do eq=1,toteq
                     call rzero(flux(1,f,e,eq),nxz)
                  enddo
               endif
            endif
         enddo
      enddo

      return
      end

      subroutine a5adiabatic_wall(eflx,f,e,dU,wstate)
      include 'SIZE'
      include 'INPUT'
      include 'GEOM' ! for UNX under ADIABATIC WALL METHOD "A"
      include 'CMTDATA'
      real eflx  (nx1*nz1,2*ndim,nelt) ! better be zero on entry
      real dU    (nx1*nz1,2*ndim,nelt,toteq)
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      common /scrns/ flxscr(lx1*lz1)
      real flxscr

      nxz=nx1*nz1

      call rzero(eflx(1,f,e),nxz)
      call rzero(hface,nxz)

      call a51dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,unx(1,1,f,e),nxz)
      call a52dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,uny(1,1,f,e),nxz)
      if (if3d) then
         call a53dUadia(flxscr,f,e,dU,wstate)
         call add2col2(eflx(1,f,e),flxsxcr,unz(1,1,f,e),nxz)
      endif
      return
      end

      subroutine a51dUadia(flux,f,ie,dU,wstate)
! same as A51 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5x+cv*lambda*u1*dU4z-kmcvmu*u3*dU4x+cv*lambda*u1*dU3y
     1   -kmcvmu*u2*dU3x+cv*mu*u3*dU2z+cv*mu*u2*dU2y+(cv*lambda-
     2   K+2*cv*mu)*u1*dU2x-cv*lambdamu*u1*u3*dU1z-cv*lambdamu
     3   *u1*u2*dU1y+(K*u3**2-cv*mu*u3**2+K*u2**2-cv*mu*u2**2-cv*la
     4   mbda*u1**2+K*u1**2-2*cv*mu*u1**2-E*K)*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a52dUadia(flux,f,ie,dU,wstate)
! same as A52 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5y+cv*lambda*u2*dU4z-kmcvmu*u3*dU4y+cv*mu*u3*dU3z+(cv
     1   *lambda-K+2*cv*mu)*u2*dU3y+cv*mu*u1*dU3x-kmcvmu*u1*dU2y+
     2   cv*lambda*u2*dU2x-cv*lambdamu*u2*u3*dU1z+(K*u3**2-cv*mu
     3   *u3**2-cv*lambda*u2**2+K*u2**2-2*cv*mu*u2**2+K*u1**2-cv*mu*
     4   u1**2-E*K)*dU1y-cv*lambdamu*u1*u2*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a53dUadia(flux,f,ie,dU,wstate)
! same as A53 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*(dU5z-E*dU1z)+c_v*u3*(lambda*dU4z+2*mu*dU4z+lambda*dU3y+lambda
     1   *dU2x)-K*u3*dU4z+c_v*mu*u2*(dU4y+dU3z)+c_v*mu*u1*(dU4x+dU2z)-
     2   K*u2*dU3z-K*u1*dU2z-c_v*(lambda+2*mu)*u3**2*dU1z+K*u3**2*dU1z+
     3   K*u2**2*dU1z-c_v*mu*u2**2*dU1z+K*u1**2*dU1z-c_v*mu*u1**2*dU1z-c
     4   _v*(lambda+mu)*u2*u3*dU1y-c_v*(lambda+mu)*u1*u3*dU1x)/(c_v*rho)
      enddo
      return
      end
