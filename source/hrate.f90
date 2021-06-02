module hrate

  implicit none

contains

  !=====================================================================
  !> @brief High level wrapper to apply cooling
  !> @details High level wrapper to apply cooling
  !! @n  parametrized cooling curve, uses the ionization state of
  !! hydrogen and ties the O I and II to it
  subroutine update_neutral_fraction()

    use parameters!, only : neq, nx, ny, nz, tsc, dif_rad, charge_exchange
    use globals!, only : u, primit, coords, dt_CFL
    use difrad!, only : ph

    implicit none
    real    :: dt_seconds
    integer :: i,j,k
    integer :: mark, nb, bID

    dt_seconds = dt*t_sc

    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        if (dif_rad) then
          call solve_h_rate(bID, dt_seconds, ph(bID,i,j,k) )
        else
          call solve_h_rate(bID, dt_seconds, 0.0)
        end if

      end if
    end do

  end subroutine update_neutral_fraction

  !======================================================================
  !> @brief calculates the recombination rate (case B)
  !> @details calculates the recombination rate (case B)
  !> @param real8 [in] T : Temperature K
  function alpha(T)

    implicit none

    real (kind=8) :: alpha
    real (kind=8), intent(in) :: T

    alpha=2.55d-13*(1.d4/T)**0.79

  end function alpha

  !======================================================================
  !> @brief calculates the collisional ionization rate
  !> @details calculates the collisional ionization rate
  !> @param real8[in] T : Temperature K
  function colf(T)

    implicit none

    real (kind=8) :: colf
    real (kind=8), intent(in) :: T

    colf=5.83d-11*sqrt(T)*exp(-157828./T)

  end function colf


  !=======================================================================
  !> @brief Updates the ionization fraction using Hrate eqn.
  !> @param real [in] dt        : timestep (seconds)
  !> @param real [in] uu(neq)   : conserved variables in one cell
  !> @param real [in] prim(neq) : primitive variables in one cell
  !> @param real [in] tau       : optical depth (not in use)
  !> @param real [in] radphi    : photoionizing rate
  subroutine solve_h_rate(bIndx,dt,radphi)

    use parameters
    use hydro_core
    implicit none

    real, intent(in)      :: dt
    integer, intent(in)   :: bIndx
    real, intent(in)      :: radphi(bIndx,:,:,:)
  !  real, intent(inout),dimension(neq) :: uu, prim
    real                               :: T
    real (kind=8) :: dh, y0, g0, e, y1
    real (kind=8) :: fpn
    real (kind=8) :: col,rec,a,b,c,d
    !      xi - neutral carbon abundance (for non-zero electron density)
    real (kind=8), parameter ::  xi=1.e-4

    do i=1,ncells_x
      do j=1,ncells_y
        do k=1,ncells_z

          !   solve for the ionization fraction and the internal energy

          call calcTemp (PRIM(bIndx,:,i,j,k), T)
          col=colf(real(T,8))               !# collisional ionization rate
          rec=alpha(real(T,8))              !# rad. recombination rate

          y0=real( U(bIndx,neqdyn+1,i,j,k)/ U(bIndx,1,i,j,k), 8 )  !# neutral H fraction
          dh=real( U(bIndx,1,i,j,k), 8 )               !# H density
          fpn=real(radphi(Bindx,i,j,k), 8)/dh            !# ionizing flux per nucleus

          !    solve for the new neutral fraction using the analytical
          !    solution (see notes)
          a=rec+col

          if (dif_rad) then
            b=-((2.+xi)*rec+(1.+xi)*col+fpn)
          else
            b=-((2.+xi)*rec+(1.+xi)*col    )
          end if

          c=(1.+xi)*rec
          d=sqrt(b**2-4.*a*c)
          g0=(2.*a*y0+b+d)/(2.*a*y0+b-d)
          e=exp( -d*dh*real(dt,8) )

          y1=(-b-d*(1.+g0*e)/(1.-g0*e))/(2.*a) !# the new neutral fraction
          y1=min(y1,1.-xi)
          y1=max(y1,0.)

          !   update the density of neutrals (conserved vars)
          U(bIndx,neqdyn+1,i,j,k)=real(y1)*U(bIndx,1,i,j,k)

        end do
      end do
    end do
  end subroutine solve_h_rate

  !======================================================================

end module hrate
