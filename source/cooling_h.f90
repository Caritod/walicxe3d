!> @brief Cooling with parametrized cooling and H rate equation
!> @details Cooling with parametrized cooling and H rate equation

module cooling_H

!#ifdef PASSIVES

  implicit none

contains

  !=======================================================================
  !> @brief High level wrapper to apply cooling
  !> @details High level wrapper to apply cooling
  !> @n  parametrized cooling curve, uses the ionization state of
  !> hydrogen and ties the O I and II to it
  subroutine coolingh()

    use parameters
    use globals
    use tictoc
    use difrad
    implicit none

    integer :: mark, nb, bID
    real :: maxloss

    if (cooling_type.ne.COOL_NONE) then

      if (verbosity > 3) call tic(mark)
      if (verbosity > 1) then
        write(logu,*) ""
        write(logu,'(1x,a)') "============================================"
        write(logu,'(1x,a)') " Applying Radiative H Cooling ..."
        write(logu,'(1x,a)') "============================================"
        write(logu,*) ""
      end if
      select case (cooling_type)

      case (COOL_H)

        ! Apply tabulated cooling to all local blocks
        do nb=1,nbMaxProc
          bID = localBlocks(nb)
          if (bID.ne.-1) then

            !call apply_cooling (nb, maxloss)
            if (dif_rad) then
              call cooling_h_neq(nb, phi(nb,:,:,:),maxloss)
              !call cooling_h_neq(primit(:,i,j,k), u(:,i,j,k), dt_seconds,ph(i,j,k))
            else
              call cooling_h_neq(nb,0,maxloss)
            end if

            if (maxloss.ge.cooling_limit) then
              write(logu,'(1x,a,i0,a,f6.3,a,f6.3)') &
                "Cooling warning for block ", bID, "! Max loss ", &
                maxloss, ", limited to", cooling_limit
            end if

          end if
        end do

      end select

      if (verbosity > 3) write(logu,'(1x,a,a)') "> Cooling applied in ", nicetoc(mark)

    end if


  end subroutine coolingh

  !======================================================================
  !> @brief betaH(T)
  !> @details @f$ \beta_H(T) @f$
  !> @param real 8[in] T : Temperature K
  function betah(T)

    implicit none

    real (kind=8) ::  betah
    real (kind=8), intent(in) ::  T
    real (kind=8)             ::  a

    a=157890./T
    betah=1.133D-24/sqrt(a)*(-0.0713+0.5*log(a)+0.640*a**(-0.33333))

  end function betah

  !======================================================================
  !> @brief Non equilibrium cooling
  !> @details   Non-equilibrium energy loss for low temperatures
  !>     considering the collisional excitation of [O I] and
  !>   [O II] lines and radiative recombination of H. This
  !>   cooling rate is multiplied by a factor of 7.033 so
  !>   that it has the same value as the "coronal equilibrium"
  !>   cooling rate at a temperature of 44770 K (at temperatures
  !>   higher than this value, the equilibrium cooling rate is
  !>   used). The collisional ionization of H and excitation
  !>   of Lyman-alpha are computed separately, and added to
  !>   the cooling rate.
  !> @param real8 [in] x1  : initial H ionization fraction
  !> @param real8 [in] x2  : final H ionization fraction
  !> @param real [in] dt  : timestep
  !> @param real8 [in] den : total density of hydrogen
  !> @param real8 [in] dh0 : density of neutral hydrogen
  !> @param real8 [in] Te0 : Temperature
  function aloss(x1,x2,den,dH0,Te0)

    implicit none

    real (kind=8)             :: aloss
!    real, intent(in)          :: dt
    real (kind=8), intent(in) :: x1,x2,den,dH0,Te0
    real, parameter :: XION=2.179e-11, XH=0.9,XO=1.e-3
    real, parameter :: C0=0.5732,C1=1.8288e-5,C2=-1.15822e-10,C3=9.4288e-16
    real, parameter :: D0=0.5856,D1=1.55083e-5,D2=-9.669e-12, D3=5.716e-19
    real, parameter :: ENK=118409.,EN=1.634E-11
    real (kind=8)   :: Te, dHp,de,dOI,dOII,omega,omegaL,omegaH,frac,qla
    real (kind=8)   :: ecoll,cion,eion,erec,Tm,T2,eOI,eOII,equil,fr,ex2,tanh
    real (kind=8)   :: betaf, HIIcool

    Te  = max(Te0,10.)
    dHp = (1.-X1)*den    ! hydrogen density
    de  = dHp+1.E-4*den  ! electron density
    dOI = XO*dH0         ! oxigen density
    dOII= XO*dHp         ! ionized oxigen density

    if(Te <= 1e4 ) then
      aloss = 1e-30
      return
    end if

    !   Collisionally excited Lyman alpha
    if(Te .le. 55000.) omega = C0 + Te*(C1 + Te*(C2 + Te*C3))
    if(Te .ge. 72000.) omega = D0 + Te*(D1 + Te*(D2 + Te*D3))
    if(Te .gt. 55000. .and. Te .lt. 72000.) then
      omegaL = C0 + Te*(C1 + Te*(C2 + Te*C3))
      omegaH = D0 + Te*(D1 + Te*(D2 + Te*D3))
      frac   = (Te-55000.)/17000.
      omega  = (1.-frac)*omegaL + frac*omegaH
    end if
    qla   = 8.6287E-6/(2.*sqrt(Te))*omega*exp(-ENK/Te)
    ecoll = de*dH0*qla*EN
    ecoll = max(ecoll,0.)

    !   Hydrogen recombination and collisional ionization
    cion = 5.834E-11*sqrt(Te)*exp(-1.579E5/Te)
    eion = de*dH0*cion*XION
    erec = de*dHp*(betah(Te))
    erec = max(erec,0.)

    !   [O I] and [O II] coll. excited lines
    Tm   = 1./Te
    Tm2  = Tm*Tm
    eOI  = de*dOI*10.**(1381465*Tm2-12328.69*Tm-19.82621)
    eOII = de*dOII*10.**(-2061075.*Tm2-14596.24*Tm-19.01402)
    eOI  = max(eOI,0.)
    eOII = max(eOII,0.)

    !   free-free cooling
    betaf  = 1.3*1.42E-27*Te**0.5
    HIIcool= de*dHp*betaf

    !   Equilibrium cooling (for high Te)
    equil = (1.0455E-18/Te**0.63)*(1. - exp(-(Te*1.E-5)**1.63))*de*den + HIIcool

    !   switch between non-equil. and equil. ionization cooling
    if(Te .le. 44770.) fr = 0.
    if(Te .ge. 54770.) fr = 1.
    if(Te .gt. 44770. .and. Te .lt. 54770.) then
      ex2  = exp(-2.*(Te-49770.)/500.)
      tanh = (1. - ex2)/(1. + ex2)
      fr   = 0.5*(1. + tanh)
    end if

    aloss = ecoll + eion + (erec + 7.033*(eOI + eOII))*(1.-fr) + equil*fr

  end function aloss

  !=======================================================================
  !> @brief
  !> @details
  !> @param real [in] uu(neq) : primitive variablas in one cell
  !> @param real [in] uu(neq) : conserved variablas in one cell
  !> @param real [in] dt      : timestep (seconds)
  !> @param real [in] radphi  : photoionizing rate
  subroutine  cooling_h_neq(bIndx, radphi, maxloss)

    use parameters   !#add energy per ionzation as a parameter
    use globals
    use hydro_core, only : calcTemp
    use constants
    use difrad

    implicit none
!    real, intent(inout) :: uu(neq), pp(neq)  ! hace falta que los defina?


    integer, intent(in):: bIndx
    real, intent(in)   :: radphi(bIndx,:,:,:)
    real, intent(out)  :: maxloss
    real, parameter    :: T_floor = 10.0

    real(kind = 8)     :: y0, y1, dh, dh0, gain, tprime, al, ce, temp, new_temp, logT
    real               ::  vel2, ETH, EK, cool_factor
    real               :: frac_loss!, metal
    integer            :: i, j, k

    maxloss = 0.0

    do i=1,ncells_x
      do j=1,ncells_y
        do k=1,ncells_z

            y0 =  real( PRIM(bIndx,neqdyn+1,i,j,k)/PRIM(bIndx,1,i,j,k), 8)  !# neutral H fraction (t0)
            y1  = real( U(bIndx,neqdyn+1,i,j,k)/U(bIndx,1,i,j,k)      , 8)  !# neutral H fraction (t0+dt) fraccion actualizada
            dh  = real( PRIM(bIndx,1,i,j,k)                           , 8)   !# total NH
            dh0 = real( PRIM(bIndx,neqdyn+1,i,j,k)                    , 8)   !# neutrals density

            ! Calculate temperature of this cell
            call calcTemp (PRIM(bIndx,:,i,j,k), temp)
            logT = log10(temp)
!              call u2prim(uu,pp,T)               !# get temperature

            ! Cooling not applied below cool_Tmin
            if (logT.ge.cool_Tmin) then

                !  get the energy losses
                al = aloss(y0,y1,dh,dh0,real(temp,8)) !chequear las unidades

                if (dif_rad) then
                  gain   = real(radphi(bIndx,i,j,k),8)*dh0*Kb*energy_per_ionization
                  tprime = max( gain*real(temp,8)/al,7000.)
                else
                  tprime=10.
                end if

                ! Calculate radiative loss and new temperature
                ! numdens = PRIM(bIndx,1,i,j,k)*d_sc/(mui*AMU)  ! cgs, gas ionized
                ! ce = aloss*(numdens**2)/(ETH*e_sc)  ! cgs
                ! cool_factor = exp(-ce*(dt*t_sc))
                ! frac_loss = 1.0-cool_factor

                ce = (2.*al/dh)/(3.*Kb*real(temp,8)) !!chequear que este en cgs
                cool_factor = exp(-ce*(dt*t_sc))   !  in physical units
                frac_loss = 1.0-cool_factor

                ! Record maximum cooling for this block before limiting
                maxloss = max(maxloss, frac_loss)

                ! Limit cool_factor directly, if needed
                if (cool_factor.lt.1.0-cooling_limit) then
                  cool_factor = 1.0-cooling_limit
                end if

                ! Impose a temperature floor by adjusting cool_factor, if needed
                ! Set to 10 K by default
                new_temp = temp * cool_factor

                if (new_temp.lt.T_floor) then
                  new_temp = T_floor
                  cool_factor = T_floor / temp
                end if

                ! Update pressure and total energy
                PRIM(bIndx,5,i,j,k) = PRIM(bIndx,5,i,j,k) * cool_factor

                ETH = CV * PRIM(bIndx,5,i,j,k)
                vel2 = PRIM(bIndx,2,i,j,k)**2 + &
                       PRIM(bIndx,3,i,j,k)**2 + &
                       PRIM(bIndx,4,i,j,k)**2
                EK = 0.5 * PRIM(bIndx,1,i,j,k) * vel2
                U(bIndx,5,i,j,k) = EK + ETH

           end if

        end do
      end do
    end do
end subroutine cooling_h_neq


!======================================================================

!#endif

end module cooling_H
