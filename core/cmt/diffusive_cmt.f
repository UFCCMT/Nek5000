      subroutine imqqtu(ummcu,uminus,uplus)
! Computes (I-0.5*QQT)U for all five conserved variables.
! See call in compute_rhs_and_dt for important documentation
!                                     -
! Spoiler: for SEM this comes out to U -{{U}}, which when
!          spoken is "UMinus Minus the Central flux of U" which I
!          then abbreviate as ummcu
      include 'SIZE'

      real ummcu (nx1*nz1*2*ndim*nelt,toteq) ! intent(out)
      real uminus(nx1*nz1*2*ndim*nelt,toteq) ! intent(in)
      real uplus (nx1*nz1*2*ndim*nelt,toteq) ! intent(in)
      integer ivar

      nf = nx1*nz1*2*ndim*nelt
      const=-0.5
! U-{{U}} on interior faces. first just do it on all faces.
      do ivar=1,toteq
         call add3(ummcu(1,ivar),uminus(1,ivar),uplus(1,ivar),nf)
         call cmult(ummcu(1,ivar),const,nf)        !         -
         call add2(ummcu(1,ivar),uminus(1,ivar),nf)!ummcu = U -{{U}}
      enddo

! v+ undefined on boundary faces, so (I-0.5QQ^T) degenerates to 
! [[U]] with respect to the Dirichlet boundary state
      call imqqtu_dirichlet(ummcu,uminus,uplus)

      return
      end

!-----------------------------------------------------------------------

      subroutine imqqtu_dirichlet(ummcu,uminus,uplus)
      return
      end

!-----------------------------------------------------------------------

      subroutine fluxj_ns(flux,gradu,e,eq)
! viscous flux jacobian for compressible Navier-Stokes equations (NS)
! SOLN and CMTDATA are indexed, assuming vdiff has been filled by uservp
! somehow. In serious need of debugging and replacement.
      include 'SIZE'
      include 'INPUT'! TRIAGE?
      include 'SOLN' ! TRIAGE?

      parameter (ldd=lx1*ly1*lz1)
      common /ctmp1/ viscscr(lx1,ly1,lz1)
      real viscscr

      integer e,eq
      real flux(nx1*ny1*nz1,ndim),gradu(nx1*ny1*nz1,toteq,ndim)
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/

      n=nx1*ny1*nz1
      call rzero(flux,n)

! This is a disaster that I might want to program less cleverly
      if (eq .lt. toteq) then ! TRIAGE. CAN'T GET AGRADU_NS to WORK
                              ! for ENERGY EQUATION. MAXIMA ROUTINES
                              ! BELOW
!!      do j=1,ndim
!!         do k=1,ndim
!!            ieijk=0
!!!           if (eq .lt. toteq .and. eq .gt. 1) ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
!!            if (eq.gt.1)ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
!!
!!            if (ieijk .eq. 0) then
!!              call agradu_ns(flux(1,j),gradu(1,1,k),viscscr,e,
!!     >                           eq,j,k)
!!            endif
!!         enddo
!!      enddo
! JH110716 Maxima routines added for every viscous flux.
!          agradu_ns has failed all verification checks for homentropic vortex
!          initialization.
!          start over
        if (eq.eq.2) then
           call A21kldUldxk(flux(1,1),gradu,e)
           call A22kldUldxk(flux(1,2),gradu,e)
           call A23kldUldxk(flux(1,3),gradu,e)
        elseif (eq.eq.3) then
           call A31kldUldxk(flux(1,1),gradu,e)
           call A32kldUldxk(flux(1,2),gradu,e)
           call A33kldUldxk(flux(1,3),gradu,e)
        elseif (eq.eq.4) then
           call A41kldUldxk(flux(1,1),gradu,e)
           call A42kldUldxk(flux(1,2),gradu,e)
           call A43kldUldxk(flux(1,3),gradu,e)
        endif

      else ! Energy equation courtesy of thoroughly-checked maxima
           ! until I can get agradu_ns working correctly
         if (if3d) then
            call a53kldUldxk(flux(1,3),gradu,e)
         else
            call rzero(gradu(1,1,3),nx1*ny1*nz1*toteq)
            call rzero(vz(1,1,1,e),nx1*ny1*nz1)
         endif
         call a51kldUldxk(flux(1,1),gradu,e)
         call a52kldUldxk(flux(1,2),gradu,e)
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine fluxj(flux,gradu,e,eq)
! JH110216 
! I might end up using this, but just for \nabla (log rho) stuff due to
! entropy viscosity. Equations (5.13-14) from Guermond & Popov (2014)
! SIAM J. Appl. Math. 74 imply, however, that his preferred regularization
! has extra off-diagonal terms in the momentum and energy equations.
! Not sure where to put those yet.
      include 'SIZE'

      integer e,eq
      real flux(nx1*ny1*nz1,ndim),gradu(nx1*ny1*nz1,toteq,ndim)

      n=nx1*ny1*nz1
      do j=1,ndim
         call agradu(flux(1,j),gradu(1,1,j),e,eq,j)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine agradu_sfc(gijklu,gvar,du,eq)
! JH110116 DEPRECATED. Always apply A to volume, not surface points.
!                      igtu_cmt will reflect this change in philosophy.
!                      Less to debug that way
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
! subroutine for computing flux of a conserved variable by higher-order
! differential operators.
! This one is classic Navier-Stokes, that is, it computes viscous fluxes
! for everything except gas density.
! eq         index i; LHS equation
! jflux      index j; flux direction
! kdir       index k; direction of derivative or jump in U
      parameter (lxyz=lx1*lz1*2*ldim)
      integer  eq,jflux
      real    du(nx1*nz1*2*ndim*nelt,toteq)
      real gvar(nx1*nz1*2*ndim*nelt,*)    ! intent(in)
! variables making up Gjkil terms, viscous stress tensor and total energy
! equation, compressible Navier-Stokes equations
! assume the following ordering remains in CMTDATA
!     gvar(:,1)  rho ! especially here
!     gvar(:,2)  u   ! especially here
!     gvar(:,3)  v   ! especially here
!     gvar(:,4)  w   ! especially here
!     gvar(:,5)  p
!     gvar(:,6)  T
!     gvar(:,7)  a
!     gvar(:,8)  phi_g
!     gvar(:,9)  rho*cv
!     gvar(:,10) rho*cp
!     gvar(:,11) mu
!     du(1,:)  rho
!     du(2,:)  rho u
!     du(3,:)  rho v
!     du(4,:)  rho w
!     du(5,:)  rho E

      real gijklu(nx1*nz1*2*ndim*nelt)

      npt=lxyz*nelt ! lazy
      call col3(gijklu,du(1,eq),gvar(1,imuf),npt)

      return
      end

!-----------------------------------------------------------------------

      subroutine agradu(gijklu,dut,e,eq,jflux)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
! monolithic viscous flux jacobian
! eq         index i; LHS equation and dut variable
! jflux      index j; flux direction
      parameter (lxyz=lx1*ly1*lz1)
      integer  e,eq,jflux
      real   dut(lxyz,toteq)
! derivatives of conserved variables gradu
!     dut(:,1) rho
!     dut(:,2) rho u
!     dut(:,3) rho v
!     dut(:,4) rho w
!     dut(:,5) rho E

      real gijklu(lxyz)

      call col3(gijklu,dut(1,eq),vdiff(1,1,1,e,imu),lxyz)

      return
      end

!-----------------------------------------------------------------------

      subroutine agradu_ns_sfc(gijklu,gvar,du,visco,eq,jflux,kdir)
! JH110116 DEPRECATED. Always apply A to volume, not surface points.
!                      igtu_cmt will reflect this change in philosophy.
!                      Less to debug that way
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
! subroutine for computing flux of a conserved variable by higher-order
! differential operators.
! This one is classic Navier-Stokes, that is, it computes viscous fluxes
! for everything except gas density.
! eq         index i; LHS equation
! jflux      index j; flux direction
! kdir       index k; direction of derivative or jump in U
      parameter (lxyz=lx1*lz1*2*ldim)
      integer  eq,jflux,kdir
      real    du(lxyz*nelt,toteq)
      real visco(lxyz*nelt) ! you know, you should probably just
                         ! pass mu+lambda and mu-k/cv when eq=5
                         ! so you don't have to recompute them
                         ! so many times
      real gvar(lxyz*nelt,*)    ! intent(in)
! variables making up Gjkil terms, viscous stress tensor and total energy
! equation, compressible Navier-Stokes equations
! assume the following ordering remains in CMTDATA
!     gvar(:,1)  rho ! especially here
!     gvar(:,2)  u   ! especially here
!     gvar(:,3)  v   ! especially here
!     gvar(:,4)  w   ! especially here
!     gvar(:,5)  p
!     gvar(:,6)  T
!     gvar(:,7)  a
!     gvar(:,8)  phi_g
!     gvar(:,9)  rho*cv
!     gvar(:,10) rho*cp
!     gvar(:,11) mu
!     gvar(:,12) thermal conductivity
!     gvar(:,13) lambda
!     gvar(:,18) U5 ! FIX THE DAMN ENERGY COMPUTATION I'M SIGHING ABOUT
! derivatives or jumps, conserved variables, compressible Navier-Stokes
! equations
!     du(1,:)  rho
!     du(2,:)  rho u
!     du(3,:)  rho v
!     du(4,:)  rho w
!     du(5,:)  rho E
      real gijklu(lxyz*nelt) !incremented. never give this exclusive intent
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/

      if (eq .eq. 1) return

      npt=lxyz*nelt ! lazy

      call rzero(visco,npt)

      if (eq .lt. 5) then

         if (jflux .eq. eq-1) then
            m=kdir
            call copy(visco,gvar(1,ilamf),npt)
            if (kdir .eq. jflux) then
               call add2s2(visco,gvar(1,imuf),2.0,npt)
            endif
         else
            call copy(visco,gvar(1,imuf),npt)
            if (kdir .eq. jflux) then
               m=eq-1
            else
               m=jflux
            endif
         endif

         m=m+1 ! skip density

         call invcol2(visco,gvar(1,irho),npt)
         call subcol4(gijklu,visco,gvar(1,m),du(1,1),npt)
         call addcol3(gijklu,visco,du(1,m),npt)

      else ! energy equation is very different. and could use a rewrite

         if (jflux .eq. kdir) then
            kp1=kdir+1

            l=1 ! sigh
            call vdot3(visco,gvar(1,iux),gvar(1,iux), ! now twoke
     >                       gvar(1,iuy),gvar(1,iuy),
     >                       gvar(1,iuz),gvar(1,iuz),npt)
            do ipt=1,npt ! someone else can make this loop more clever
               gdu=(gvar(ipt,imuf)-gvar(ipt,ikndf))*twoke
               energy=gvar(ipt,icvf)*gvar(ipt,ithm)+0.5*twoke ! sigh. iu5?/phi?
               gdu=gdu+(gvar(ipt,imuf)+gvar(ipt,ilamf))*gvar(ipt,kp1)**2
               gijklu(ipt)=gijklu(ipt)-(gvar(ipt,ikndf)*energy-gdu)*
     >                                  du(ipt,l)/gvar(ipt,irho)
            enddo

            call sub3(visco,gvar(1,imuf),gvar(1,ikndf),npt) ! form mu-K/cv
            do ipt=1,npt
               visco(ipt)=visco(ipt)/gvar(ipt,1)
               gdu=0.0
               do l=2,ldim+1 ! both gvar and du are indexed by l
                  gdu=gdu+gvar(ipt,l)*du(ipt,l)
               enddo
               gijklu(ipt)=gijklu(ipt)+gdu*visco(ipt)
            enddo
            call add3(visco,gvar(1,imuf),gvar(1,ilamf),npt)
            l=jflux+1
            do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)+visco(ipt)*gvar(ipt,l)*du(ipt,l)/
     >                                                      gvar(ipt,1)
            enddo

            l=5
            call copy(visco,gvar(1,ikndf),npt)
            call invcol2(visco,gvar(1,1),npt)
            call add2col2(gijklu,visco,du(1,l),npt)

         else ! dU is off-diagonal

            call add3(visco,gvar(1,imuf),gvar(1,ilamf),npt)
            jp1=jflux+1
            kp1=kdir+1
            call invcol2(visco,gvar(1,1),npt)
            do ipt=1,npt
               gijklu(ipt)=gijklu(ipt)-
     >                  visco(ipt)*gvar(ipt,jp1)*gvar(ipt,kp1)*du(ipt,1)
            enddo

            do l=2,ldim+1
               lm1=l-1
               if (eijk3(jflux,kdir,lm1) .eq. 0) then
                  if (lm1 .eq. kdir) then
                     call copy(visco,gvar(1,ilamf),npt)
                     m=jflux+1
                  else
                     call copy(visco,gvar(1,imuf),npt)
                     m=kdir+1
                  endif
                  call invcol2(visco,gvar(1,1),npt)
                  call addcol4(gijklu,visco,gvar(1,m),du(1,l),npt)
               endif
            enddo ! l

         endif ! diagonal?

      endif ! energy equation

      return
      end

!-----------------------------------------------------------------------

      subroutine agradu_ns(gijklu,dut,visco,e,eq,jflux,kdir)
! JH110716. Declaring defeat. Discard and start over
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'CMTDATA'
      include 'GEOM' ! diagnostic
! classic Navier-Stokes flux jacobian that computes viscous fluxes
! for everything except gas density.
! eq         index i; LHS equation
! jflux      index j; flux direction
! kdir       index k; direction of derivative or jump in U
      parameter (lxyz=lx1*ly1*lz1)
      integer  e,eq,jflux,kdir
      real   dut(lxyz,toteq)
      real visco(lxyz) ! you know, you should probably just
                         ! pass mu+lambda and mu-k/cv when eq=5
                         ! so you don't have to recompute them
                         ! so many times
! derivatives or jumps, conserved variables, compressible Navier-Stokes
! equations
!     du(1,:) or dut(:,1) rho
!     du(2,:) or dut(:,2) rho u
!     du(3,:) or dut(:,3) rho v
!     du(4,:) or dut(:,4) rho w
!     du(5,:) or dut(:,5) rho E
      real gijklu(lxyz) !incremented. never give this exclusive intent
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/
      real mu,lambda,kond,cv,uiui(lxyz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /CTMP0/ and /CTMP1/ IN USE BY CALLING SUBROUTINE; NO TOUCHY!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (eq .eq. 1) return ! do this in plain old agradu

      npt=lxyz !%s/npt/lxyz/g is hard

      call rzero(visco,npt)
      call rzero(uiui,npt)

      if (eq .lt. 5) then ! momentum fluxes via factorized viscous
                          ! stress tensor

         if (jflux .eq. eq-1) then
            m=kdir
            call copy(visco,vdiff(1,1,1,e,ilam),npt)
            if (kdir .eq. jflux) then
               call add2s2(visco,vdiff(1,1,1,e,imu),2.0,npt)
            endif
         else
            call copy(visco,vdiff(1,1,1,e,imu),npt)
            if (kdir .eq. jflux) then
               m=eq-1
            else
               m=jflux
            endif
         endif

         m=m+1

         call invcol2(visco,vtrans(1,1,1,e,irho),npt)
         call addcol3(gijklu,visco,dut(1,m),npt)
! sigh. times like this I hate fortran. AdU-= drho/dx_k contrib
         if (m .eq. 2) call subcol4(gijklu,visco,dut,vx(1,1,1,e),npt)
         if (m .eq. 3) call subcol4(gijklu,visco,dut,vy(1,1,1,e),npt)
         if (m .eq. 4) call subcol4(gijklu,visco,dut,vz(1,1,1,e),npt)

      else ! energy equation is very different. and could use a rewrite

         if (jflux .eq. kdir) then
! first dU1 flux, then dU2-dU3 for work term, then dU5
            if (if3d) then ! uiui= 2KE/rho=u_i u_i
               call vdot3(uiui,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),
     >                         vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),npt)
            else
               call vdot2(uiui,vx(1,1,1,e),vy(1,1,1,e),
     >                         vx(1,1,1,e),vy(1,1,1,e),npt)
            endif
            call invcol3(visco,vdiff(1,1,1,e,imu),
     >                        vtrans(1,1,1,e,irho),npt)
            call chsign(visco,npt)
            do i=1,npt
               visco(i)=visco(i)+0.5*vdiff(i,1,1,e,iknd)/
     >                              vtrans(i,1,1,e,icv) !rho*cv
            enddo

            call addcol3(gijklu,visco,dut,npt) !*dU1, but one more to go

            call add3(visco,vdiff(1,1,1,e,imu),vdiff(1,1,1,e,ilam),npt)
            call invcol2(visco,vtrans(1,1,1,e,irho),npt)
            if(jflux .eq. 1) call col3(uiui,vx(1,1,1,e),vx(1,1,1,e),npt)
            if(jflux .eq. 2) call col3(uiui,vy(1,1,1,e),vy(1,1,1,e),npt)
            if(jflux .eq. 3) call col3(uiui,vz(1,1,1,e),vz(1,1,1,e),npt)
            do i=1,npt ! someone else can vectorize this better
               temp=t(i,1,1,e,1)
               rho=vtrans(i,1,1,e,irho)
               kond=vdiff(i,1,1,e,iknd)
               gijklu(i)=gijklu(i)-(kond*temp/rho+visco(i)*uiui(i))*
     >                   dut(i,1)
            enddo ! done with *dU1 flux

            call invcol3(visco,vdiff(1,1,1,e,imu), ! visco=mu/rho
     >                       vtrans(1,1,1,e,irho),npt)
            call invcol3(uiui,vdiff(1,1,1,e,iknd), ! uiui=krcv
     >                        vtrans(1,1,1,e,icv),npt) ! rho*cv
            call sub2(visco,uiui,npt) ! visco=1/rho(mu-k/cv)
            call addcol4(gijklu,visco,dut(1,2),vx(1,1,1,e),npt)
            call addcol4(gijklu,visco,dut(1,3),vy(1,1,1,e),npt)
            if (if3d)call addcol4(gijklu,visco,dut(1,4),vz(1,1,1,e),npt)

            call add3(visco,vdiff(1,1,1,e,imu),vdiff(1,1,1,e,ilam),npt)
            m=jflux+1
           if(m.eq.2)call addcol4(gijklu,visco,dut(1,m),vx(1,1,1,e),npt)
           if(m.eq.3)call addcol4(gijklu,visco,dut(1,m),vy(1,1,1,e),npt)
           if(m.eq.4)call addcol4(gijklu,visco,dut(1,m),vz(1,1,1,e),npt)

            call addcol3(gijklu,uiui,dut(1,toteq),npt) !krcv*dU5

         else ! dUeverythingelse
! NOTATION ABUSE!!!!! uiui=u_{jflux}, visco=u_{kdir} !!!!!!!!!!
! math.f not to be used for actual functions of transport coefs
            if (jflux.eq.1) call copy(uiui, vx(1,1,1,e),npt)
            if (jflux.eq.2) call copy(uiui, vy(1,1,1,e),npt)
            if (jflux.eq.3) call copy(uiui, vz(1,1,1,e),npt)
            if (kdir .eq.1) call copy(visco,vx(1,1,1,e),npt) 
            if (kdir .eq.2) call copy(visco,vy(1,1,1,e),npt) 
            if (kdir .eq.3) call copy(visco,vz(1,1,1,e),npt) 
            do i=1,npt ! again, someone else can vectorize this
               mu    =vdiff (i,1,1,e,imu)
               lambda=vdiff (i,1,1,e,ilam)
               rrho  =vtrans(i,1,1,e,irho)
               gijklu(i)=gijklu(i)-(mu+lambda)*rrho*dut(i,1)*
     >                                         uiui(i)*visco(i) ! ujuk
               gijklu(i)=gijklu(i)+mu*visco(i)*rrho*dut(i,jflux+1)
               gijklu(i)=gijklu(i)+lambda*uiui(i)*rrho*dut(i,kdir+1)
            enddo
         endif

      endif ! energy equation

      return
      end

!-----------------------------------------------------------------------

      subroutine half_iku_cmt(res,diffh,e)
      include 'SIZE'
      include 'MASS'
! diffh has D AgradU. half_iku_cmt applies D^T BM1 to it and increments
! the residual res with the result
      integer e ! lopsided. routine for one element must reference bm1
      real res(nx1,ny1,nz1),diffh(nx1*ny1*nz1,ndim)

      n=nx1*ny1*nz1

      do j=1,ndim
         call col2(diffh(1,j),bm1(1,1,1,e),n)
      enddo

!     const=-1.0 ! I0
      const=1.0  ! *-1 in time march
      call gradm11_t(res,diffh,const,e)

      return
      end

!-----------------------------------------------------------------------

      subroutine compute_transport_props
! get vdiff props (viscosity in imu, second viscosity in ilam, and
! thermal conductivity in iknd; second viscosity is usually -2/3
! visc, but we refuse to assume Stokes' hypothesis for the user)
! via nekasn
! JH082216 Guermond eddy viscosity method (EVM) regularization starts
!          here
      include 'SIZE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'SOLN'
      include 'CMTDATA'

      integer   e

      do e=1,nelt
         ieg=lglel(e)
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            call nekasgn(i,j,k,e)
            call uservp(i,j,k,ieg)
            vdiff(i,j,k,e,imu) = mu
            vdiff(i,j,k,e,ilam) = lambda
            vdiff(i,j,k,e,iknd) = udiff
         enddo
         enddo
         enddo
      enddo
      return
      end

!-----------------------------------------------------------------------
! TRIAGE BELOW UNTIL I CAN FIX AGRADU_NS
!-----------------------------------------------------------------------
      subroutine a51kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
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

      subroutine a52kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
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
      subroutine a53kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
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

      subroutine A21kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=
     >(lambda*(dU4z+dU3y-u3*dU1z-u2*dU1y)+lambdamu*(dU2x-u1*dU1x))/rho
      enddo
      return
      end
      subroutine A22kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU3x=dU(i,3,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         flux(i)=mu*(dU3x+dU2y-u1*dU1y-u2*dU1x)/rho
      enddo
      return
      end
      subroutine A23kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU4x=dU(i,4,1)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4x+dU2z-u1*dU1z-u3*dU1x)/rho
      enddo
      return
      end

      subroutine A31kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU3x=dU(i,3,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         flux(i)=mu*(dU3x+dU2y-u1*dU1y-u2*dU1x)/rho
      enddo
      return
      end
      subroutine A32kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=(lambda*(dU4z+dU2x-u3*dU1z-u1*dU1x)+
     >   lambdamu*(dU3y-u2*dU1y))/rho
      enddo
      return
      end
      subroutine A33kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1y=dU(i,1,2)
         dU4y=dU(i,4,2)
         dU1z=dU(i,1,3)
         dU3z=dU(i,3,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4y+dU3z-u2*dU1z-u3*dU1y)/rho
      enddo
      return
      end

      subroutine A41kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU4x=dU(i,4,1)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4x+dU2z-u1*dU1z-u3*dU1x)/rho
      enddo
      return
      end
      subroutine A42kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1y=dU(i,1,2)
         dU4y=dU(i,4,2)
         dU1z=dU(i,1,3)
         dU3z=dU(i,3,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4y+dU3z-u2*dU1z-u3*dU1y)/rho
      enddo
      return
      end
      subroutine A43kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real lambda,mu,cv,K,rho,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=(lambda*(dU3y+dU2x-u2*dU1y-u1*dU1x)+
     >lambdamu*(dU4z-u3*dU1z))/rho
      enddo
      return
      end
