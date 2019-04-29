!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 10p3_WB_Galileon.f90
!! This file contains the definition of the WB Galileon model.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the WB Galileon model.
!! Please refer to the numerical notes for details.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_full_WB_Galileon

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use equispaced_linear_interpolation_1D
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_WB_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the WB Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_WB_Galileon

        ! the model parameters:
        real(dl)  :: c2      !< WB Galileon model parameter \f$c_2\f$
        real(dl)  :: c3      !< WB Galileon model parameter \f$c_3\f$
        real(dl)  :: c4      !< WB Galileon model parameter \f$c_4\f$
        real(dl)  :: c5      !< WB Galileon model parameter \f$c_5\f$
        real(dl)  :: XDS      !!!!!!!csi                     !< WB Galileon background parameter \f$\xi\f$ deriving from the tracker solution @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        real(dl)  :: p                 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        real(dl)  :: s                 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        real(dl)  :: Hds               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated function c (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma1      !< The interpolated function \f$\gamma_1\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma2      !< The interpolated function \f$\gamma_2\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma3      !< The interpolated function \f$\gamma_3\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma4      !< The interpolated function \f$\gamma_4\f$ (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTgamma5      !< The interpolated function \f$\gamma_4\f$ (and derivatives). @ New line @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ! some designer parameters:
        integer  :: designer_num_points = 1000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-10._dl))          !< log(a start)
        real(dl) :: x_final             = log(2._dl)                      !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBWBGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBWBGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBWBGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBWBGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBWBGalileonInitBackground               !< subroutine that initializes the background of WB Galileon.
        procedure :: solve_background                => EFTCAMBWBGalileonSolveBackground            !< subroutine that solves the background equations.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBWBGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBWBGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBWBGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBWBGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBWBGalileonParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBWBGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBWBGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBWBGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBWBGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBWBGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_WB_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBWBGalileonReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_WB_Galileon)         :: self   !< the base class
        type(TIniFile)                          :: Ini    !< Input ini file

    end subroutine EFTCAMBWBGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBWBGalileonAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_WB_Galileon) :: self !< the base class

    end subroutine EFTCAMBWBGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    !! Nothing needs to be done but procedure present because it is deferred.
    subroutine EFTCAMBWBGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_WB_Galileon)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in)     :: array  !< input array with the values of the parameters of the model.

        ! read xi:
        self%XDS = array(1) !self%csi = array(1)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%p = array(2)                      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%s = array(3)                      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    end subroutine EFTCAMBWBGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBWBGalileonInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_WB_Galileon)    :: self   !< the base class
        type(TIniFile)                     :: Ini    !< Input ini file

        self%XDS = Ini_Read_Double_File( Ini, 'WB_Galileon_XDS', 0._dl ) !self%csi = Ini_Read_Double_File( Ini, 'WB_Galileon_xi', 0._dl )@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%p = Ini_Read_Double_File( Ini, 'WB_Galileon_p', 0._dl )
        self%s = Ini_Read_Double_File( Ini, 'WB_Galileon_s', 0._dl )

    end subroutine EFTCAMBWBGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of WB Galileon.
    subroutine EFTCAMBWBGalileonInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        real(dl) :: Omega_!Omega_phi0

        Omega_phi0 = params_cache%omegav

        ! WB -> just c_4 and c_5
        self%c4 = 0.084852._dl   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%c5 = 0._dl          !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%c2 = 3._dl - 6._dl*self%c4 - 24._dl*self%c5 + 12._dl*self%c4*self%p + 24._dl*self%c4*self%q                      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        self%c3 = (sqrt(2._dl)*self%p)/(-1._dl + 2._dl*self%p + 2._dl*self%q) - (4._dl*sqrt(2._dl)*self%c4*self%p)/(-1._dl + 2._dl*self%p + 2._dl*self%q) - ( 16._dl*sqrt(2._dl)*self%WBGalileon_c5*self%p)/(-1._dl + 2._dl*self%p + 2._dl*self%q) + ( 8._dl*sqrt(2._dl)*self%c4*self%p**2._dl)/(-1._dl + 2._dl**self%p + 2._dl*self%q) -( 4._dl*sqrt(2._dl)*self%c4*self%q)/(-1._dl + 2._dl*self%p + 2._dl*self%q) - ( 16._dl*sqrt(2._dl)*self%c5*self%q)/(-1._dl + 2._dl*self%p + 2._dl*self%q) + ( 24._dl*sqrt(2._dl)*self%c4*self%p*self%q)/(-1._dl + 2._dl*self%p + 2._dl*self%q) + ( 16*sqrt(2._dl)*self%c4*self%q**2._dl)/(-1._dl + 2._dl*self%p + 2._dl*self%q)  !@@@@@@@@@@@@@@@@@@@

        self%Hds = eft_par_cache%h0_Mpc*(1-( eft_par_cache%omegac +eft_par_cache%omegab + eft_par_cache%omegag +eft_par_cache%omegar) )**(1/(2+2*self%s))

        call self%feedback()

        ! initialize interpolating functions:
        call self%EFTOmega%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTLambda%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize      ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma1%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma2%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma3%initialize ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTgamma4%initialize ( self%designer_num_points, self%x_initial, self%x_final )

        ! solve the background equations and store the solution:
        call self%solve_background( params_cache, success=success )

        success = .true.

    end subroutine EFTCAMBWBGalileonInitBackground

    subroutine EFTCAMBWBGalileonSolveBackground( self, params_cache, success )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                        :: success       !< whether the calculation ended correctly or not
        real(dl) :: PPlus, yPlus, CoeffA_Part, yStar, x
        integer  :: i
        real(dl) :: t1, t2

        t1  = self%EFTOmega%x(1)
        call output(params_cache,  1, t1 )
        do i=1, self%EFTOmega%num_points-1

            t2 = self%EFTOmega%x(i+1)
            call output(params_cache,  i+1, t2)

        end do

        success =.true.
        return

    contains

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that takes the solution of the background and stores the values of the EFT functions.
      subroutine output( eft_par_cache, ind, x)

         implicit none

         type(EFTCAMB_parameter_cache), intent(in):: eft_par_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
         integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.

         !real(dl) :: a,Omega_tot,Omega_tot_prime,Omega_tot_primeprime, Omega_tot_primeprimeprime,Omega_tot_primeprimeprimeprime
         !real(dl) :: rhonu_tot, presnu_tot, presnudot_tot, presnudotdot_tot,presnudotdotdot_tot,Omega_phi0
         !real(dl) :: rhonu, presnu, grhormass_t,presnudotdotdot,presnudotdot,presnudot, psi,psiprime, psiprimeprime, psiprimeprimeprime
         !real(dl) :: adotoa, Hdot,Hdotdot,Hddd, Hdddd, psiprimeprimeprimeprime
         !real(dl) :: phip1, phip2, phip3, phip4, phip5, m0, a2, c3
         !integer  :: nu_i

         real(dl) :: a, a2

         real(dl) :: rhonu_tot, presnu_tot, presnudot_tot, presnudotdot_tot,presnudotdotdot_tot,Omega_phi0
         real(dl) :: rhonu, presnu, grhormass_t,presnudotdotdot,presnudotdot,presnudot,
         real(dl) :: adotoa, Hdot,Hdotdot,Hddd, Hdddd

         real(dl) :: adotoaPrime, adotoaPrimePrime, adotoaPrimePrimePrime!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         real(dl) :: Phi, PhiPrime, PhiPrimePrime, PhiPrimePrimePrime, PhiPrimePrimePrimePrime, m0, a2, c3
         integer  :: nu_i

         a = Exp(x)
         a2 = a*a
         m0=1._dl

         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
           do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

             rhonu           = 0._dl
             presnu          = 0._dl
             grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

             call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
             rhonu_tot  = rhonu_tot + grhormass_t*rhonu
             presnu_tot = presnu_tot + grhormass_t*presnu

            end do
         end if

         !Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +rhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2) !@@@@@@@@@@
         adotoa =  ! solve numerical equation  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         adotoaPrime =   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         adotoaPrimePrime =   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         adotoaPrimePrimePrime =   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         Hdot = adotoaPrime*adotoa*a  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         Hdotdot = adotoaPrimePrime*adotoa**2*a**2 + adotoaPrime**2*adotoa*a**2 + adotoaPrime*adotoa**2*a  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         !Omega_phi0 = eft_par_cache%omegav
         !Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                                  & -(rhonu_tot+presnu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)
         !Hdot = adotoa**2 +0.25_dl*(eft_par_cache%h0_Mpc)**2*a**3*( 1._dl + Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_prime

         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         presnudot_tot = 0._dl
         presnudotdot_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

              rhonu           = 0._dl
              presnu          = 0._dl
              presnudot       = 0._dl
              presnudotdot    = 0._dl
              grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

              call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
              rhonu_tot  = rhonu_tot + grhormass_t*rhonu
              presnu_tot = presnu_tot + grhormass_t*presnu
              presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
              presnudotdot = eft_par_cache%Nu_pidotdot(a*eft_par_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)
              presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
              presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

            end do
         end if

         Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                                  & +(4._dl*(rhonu_tot+presnu_tot)-presnudot_tot/adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)

         !Hdotdot = 2._dl*adotoa*Hdot +3._dl*adotoa*( Hdot -adotoa**2 ) +0.25_dl*(eft_par_cache%h0_Mpc)**2*adotoa*a2**2&
                   !    & *( ( 1._dl +Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_primeprime +Omega_tot_prime**2&
                   !    & *( 4._dl*Omega_phi0/( Omega_tot**2 +4._dl*Omega_phi0 )**( 1.5_dl ) ) )

         rhonu_tot  = 0._dl
         presnu_tot = 0._dl
         presnudot_tot = 0._dl
         presnudotdot_tot = 0._dl
         presnudotdotdot_tot = 0._dl
         if ( eft_par_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, eft_par_cache%Nu_mass_eigenstates

              rhonu           = 0._dl
              presnu          = 0._dl
              presnudot       = 0._dl
              presnudotdot    = 0._dl
              presnudotdotdot = 0._dl
              grhormass_t = eft_par_cache%grhormass(nu_i)/a**2

              call eft_par_cache%Nu_background(a*eft_par_cache%nu_masses(nu_i),rhonu,presnu)
              presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
              presnudotdot = eft_par_cache%Nu_pidotdot(a*eft_par_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)
              rhonu_tot  = rhonu_tot + grhormass_t*rhonu
              presnu_tot = presnu_tot + grhormass_t*presnu
              presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
              presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))
              presnudotdotdot_tot = presnudotdotdot_tot + grhormass_t*( presnudotdotdot -12._dl*adotoa*presnudotdot &
                  & -64._dl*adotoa**3*presnu -12._dl*Hdot*presnudot +48._dl*adotoa**2*presnudot -4._dl*Hdotdot*presnu +48._dl*adotoa*Hdot*presnu)

            end do
        end if

        !Omega_tot_primeprimeprime = -60._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-6) -120._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-7)&
                    & +(-20._dl*(rhonu_tot+presnu_tot) + (6._dl/adotoa +Hdot/adotoa**3)*presnudot_tot  -1._dl/adotoa**2*presnudotdot_tot)/(eft_par_cache%h0_Mpc**2*a**5)

        Hddd = a*adotoa**3*adotoaPrime + 4*a**2*adotoa**2*adotoaPrime**2 + a**3*adotoa*adotoaPrime**3 + 3*a**2*adotoa**3*adotoaPrimePrime + 4*a**3*adotoa**2*adotoaPrime*adotoaPrimePrime + a**3*adotoa**3*adotoaPrimePrimePrime

        !Omega_tot_primeprimeprimeprime = 360._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-7) +840._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-8)&
                    & +(120._dl*(rhonu_tot+presnu_tot) + (-38._dl/adotoa -9._dl*Hdot/adotoa**3 +Hdotdot/adotoa**4 &
                    & -3._dl*Hdot**2/adotoa**5 )*presnudot_tot +presnudotdot_tot*( 9._dl/adotoa**2 +3._dl*Hdot/adotoa**4 )&
                    & -presnudotdotdot_tot/adotoa**3)/(eft_par_cache%h0_Mpc**2*a**6)

        Hdddd = a*adotoa**4*adotoaPrime + 11*a**2*adotoa**3*adotoaPrime**2 + 11*a**3*adotoa**2*adotoaPrime**3 + a**4*adotoa*adotoaPrime**4 + 7*a**2*adotoa**4*adotoaPrimePrime + 9*a**3*adotoa**3*adotoaPrime*adotoaPrimePrime + 11*a**4*adotoa**2*adotoaPrime**2*adotoaPrimePrime + 4*a**4*adotoa**3*adotoaPrimePrime**2 + 6*a**3*adotoa**4*adotoaPrimePrimePrime + 7*a**4*adotoa**3*adotoaPrime*adotoaPrimePrimePrime + a**4*adotoa**4*adotoaPrimePrimePrimePrime


          if ( a == 0._dl ) then
              return
          else if ( adotoa  == 0._dl ) then
              if  ( adotoa  == 0._dl ) return
              if  ( adotoaPrime    == 0._dl ) return !Hdot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
              if  ( adotoaPrimePrime == 0._dl ) return !Hdotdot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
              if  ( adotoaPrimePrimePrime == 0._dl ) return !Hdotdot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
              if  ( adotoaPrimePrimePrimePrime == 0._dl ) return !Hdotdot @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          end if
          !
          ! compute the psi field and its derivatives
          !psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/adotoa**2
          !psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/adotoa**4*( adotoa**2 -Hdot )
          !psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/adotoa**4*( adotoa**2 -3._dl*Hdot &
          !          & +4._dl*( Hdot/adotoa )**2 -Hdotdot/adotoa )
          !psiprimeprimeprime = -4._dl*Hdot*psiprimeprime/(a*adotoa**2)+ 2._dl*self%csi*eft_par_cache%h0_Mpc**2/(a*adotoa**5)&
          !          & *( 2._dl*adotoa*Hdot -3._dl*Hdotdot +9._dl*Hdot*Hdotdot/adotoa**2 &
          !          & -8._dl*Hdot**3/adotoa**3-Hddd/adotoa )
          !psiprimeprimeprimeprime = -4._dl/(a2*adotoa**4)*psiprimeprime*(Hdotdot*adotoa +3._dl*Hdot**2)&
          !          &- psiprimeprimeprime/a*(1._dl +9._dl*Hdot/adotoa**2) +2._dl*self%csi*eft_par_cache%h0_Mpc**2/(a2*adotoa**6)*(2._dl*Hdot**2 &
          !          &+2._dl*adotoa*Hdotdot -3._dl*Hddd+ 9._dl*Hdotdot**2/adotoa**2 -42._dl*Hdot**2*Hdotdot/adotoa**3 &
          !          &+24._dl*(Hdot/adotoa)**4 +10._dl*Hdot*Hddd/adotoa**2 -Hdddd/adotoa)
          !


          ! scalar field conversion
          !phip1 = psi/a
          !phip2 = ( psiprime*a -psi )/a2
          !phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3
          !phip4 = ( psiprimeprimeprime -3._dl*phip3 )/a
          !phip5 = -( psiprimeprimeprime -3._dl*phip3 )/a2+( psiprimeprimeprimeprime -3._dl*phip4 )/a
          !
          Chi = XDS*((a*Hds)/adotoa)**((2*s)/p)! from tracker@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         
          ChiPrime =(2*s*Chi*(adotoa - a*adotoaPrime))/(a*p*adotoa) !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          ChiPrimePrime = (2*s*Chi*(-((p - 2*s)*adotoa**2) + a**2*(p + 2*s)*adotoaPrime**2 - a*adotoa*(4*s*adotoaPrime + a*p*adotoaPrimePrime)))/(a**2*p**2*adotoa**2)!@@@@@@@@@@@@@@@@

          ChiPrimePrimePrime = (2*s*Chi*(2*(p**2 - 3*p*s + 2*s**2)*adotoa**3 - 2*a**3*(p**2 + 3*p*s + 2*s**2)*adotoaPrime**3 + 3*a**2*(p + 2*s)*adotoa*adotoaPrime*(2*s*adotoaPrime + a*p*adotoaPrimePrime) - a*adotoa**2*(-6*(p - 2*s)*s*adotoaPrime + a*p*(6*s*adotoaPrimePrime + a*p*adotoaPrimePrimePrime))))/(a**3*p**3*adotoa**3)!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          ChiPrimePrimePrimePrime = (2*s*Chi*((-6*p**3 + 22*p**2*s - 24*p*s**2 + 8*s**3)*adotoa**4 + 2*a**4*(3*p**3 + 11*p**2*s + 12*p*s**2 + 4*s**3)*adotoaPrime**4 -4*a**3*(p**2 + 3*p*s + 2*s**2)*adotoa*adotoaPrime**2*(4*s*adotoaPrime + 3*a*p*adotoaPrimePrime) + a**2*(p + 2*s)*adotoa**2*(-12*(p - 2*s)*s*adotoaPrime**2 + 3*a**2*p**2*adotoaPrimePrime**2 + 4*a*p*adotoaPrime*(6*s*adotoaPrimePrime + a*p*adotoaPrimePrimePrime)) - a*adotoa**3*(16*s*(p**2 - 3*p*s + 2*s**2)*adotoaPrime + a*p*(-12*(p - 2*s)*s*adotoaPrimePrime +  a*p*(8*s*adotoaPrimePrimePrime + a*p*adotoaPrimePrimePrimePrime)))))/(a**4*p**4*adotoa**4)!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          
          PhiPrime = Sqrt(-Chi)/adotoa    ! from Chi definition @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          PhiPrimePrime = (2*Chi*adotoaPrime - adotoa*ChiPrime)/(2.*adotoa**2*Sqrt(-Chi))           !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          PhiPrimePrimePrime =  -(adotoa**2*ChiPrime**2 + Chi**2*(-8*adotoaPrime**2 + 4*adotoa*adotoaPrimePrime) - 2*adotoa*Chi*(-2*adotoaPrime*ChiPrime + adotoa*ChiPrimePrime))/(4.*adotoa**3*(-Chi)**1.5)           !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          PhiPrimePrimePrimePrime =  (-3*adotoa**3*ChiPrime**3 + 6*adotoa**2*Chi*ChiPrime*(-(adotoaPrime*ChiPrime) + adotoa*ChiPrimePrime) + 8*Chi**3*(6*adotoaPrime**3 - 6*adotoa*adotoaPrime*adotoaPrimePrime + adotoa**2*adotoaPrimePrimePrime) - 4*adotoa*Chi**2*(6*adotoaPrime**2*ChiPrime - 3*adotoa*adotoaPrime*ChiPrimePrime + adotoa*(-3*ChiPrime*adotoaPrimePrime + adotoa*ChiPrimePrimePrime)))/(8.*adotoa**4*(-Chi)**2.5)!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! compute the background EFT functions:
          
         ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
         ! modification starts
         ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          self%EFTOmega%y(ind)    = -1 - (2*c4*adotoa**((2*p*(1 + s))/s)*PhiPrime**((2*p*(1 + s))/s))/((-1)**(p*(1 + 1/s))*XDS**(p*(1 + 1/s)))

          self%EFTOmega%yp(ind)    = (-4*c4*p*(1 + s)*adotoa**(-1 + (2*p*(1 + s))/s)*PhiPrime**(-1 + (2*p*(1 + s))/s)*(adotoa*PhiPrimePrime + PhiPrime*adotoaPrime))/((-1)**(p*(1 + 1/s))*s*XDS**(p*(1 + 1/s)))

          self%EFTOmega%ypp(ind)   = (-4*c4*p*(1 + s)*adotoa**(-2 + (2*p*(1 + s))/s)*PhiPrime**(-2 + (2*p*(1 + s))/s)*(adotoa**2*((-s + 2*p*(1 + s))*PhiPrimePrime**2 + s*PhiPrime*PhiPrimePrimePrime) + (-s + 2*p*(1 + s))*PhiPrime**2*adotoaPrime**2 + adotoa*PhiPrime*(4*p*(1 + s)*PhiPrimePrime*adotoaPrime + s*PhiPrime*adotoaPrimePrime)))/((-1)**(p*(1 + 1/s))*s**2*XDS**(p*(1 + 1/s)))

          self%EFTOmega%yppp(ind)  = (-4*c4*p*(1 + s)*adotoa**(-3 + (2*p*(1 + s))/s)*PhiPrime**(-3 + (2*p*(1 + s))/s)*(adotoa**3*(2*(s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2)*PhiPrimePrime**3 + 3*s*(-s + 2*p*(1 + s))*PhiPrime*PhiPrimePrime*PhiPrimePrimePrime + s**2*PhiPrime**2*PhiPrimePrimePrimePrime) + 2*(s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2)*PhiPrime**3*adotoaPrime**3 + 3*(-s + 2*p*(1 + s))*adotoa*PhiPrime**2*adotoaPrime*(2*p*(1 + s)*PhiPrimePrime*adotoaPrime + s*PhiPrime*adotoaPrimePrime) + adotoa**2*PhiPrime*(6*p*(1 + s)*(-s + 2*p*(1 + s))*PhiPrimePrime**2*adotoaPrime + 6*p*s*(1 + s)*PhiPrime*PhiPrimePrime*adotoaPrimePrime + s*PhiPrime*(6*p*(1 + s)*PhiPrimePrimePrime*adotoaPrime + s*PhiPrime*adotoaPrimePrimePrime))))/((-1)**(p*(1 + 1/s))*s**3*XDS**(p*(1 + 1/s)))

          self%EFTc%y(ind)         = (-2*a**2*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*p*s**2*XDS**(p + (3*p)/(2.*s))*adotoa**(2*p)*PhiPrime**(2 + 2*p) + Sqrt(2)*a*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*(3*PhiPrime - a*PhiPrimePrime) - 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*((-28*c5*s**2 + c4*p*(1 + s)*(-s + 6*p*(1 + s)))*PhiPrime**(2*(1 + p + p/s)) - a**2*c4*p*(1 + s)*(-s + 2*p*(1 + s))*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2 - a*p*(1 + s)*PhiPrime**(1 + (2*p*(1 + s))/s)*((-3*c4*s - 8*c5*s + 4*c4*p*(1 + s))*PhiPrimePrime + a*c4*s*PhiPrimePrimePrime)) - Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*adotoaPrime + 8*(-1)**p*(0,1)**(p/s)*a**2*c4*p**2*(1 + s)**2*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(2*(1 + p + p/s))*adotoaPrime**2 + 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(a*c4*p*(1 + s)*(s + 4*p*(1 + s))*PhiPrimePrime*adotoaPrime + PhiPrime*((-4*c5*s*(s + 2*p*(1 + s)) + c4*p*(1 + s)*(-s + 4*p*(1 + s)))*adotoaPrime + a*c4*p*s*(1 + s)*adotoaPrimePrime)))/(2.*(-1)**(p*(2 + 1/s))*(0,1)**(p/s)*s**2*XDS**(p*(2 + 3/(2.*s)))*PhiPrime**2)

          self%EFTLambda%y(ind)    =         ((-1)**(1 + p + p/s)*(0,1)**(p/s)*a**2*c2*Hds**2*s**2*XDS**(p + (3*p)/(2.*s))*adotoa**(2*p)*PhiPrime**(2 + 2*p) - Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*PhiPrimePrime + 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*(s*(-2*c5*s + c4*p*(1 + s))*PhiPrime**(2*(1 + p + p/s)) + a**2*c4*p*(1 + s)*(-s + 2*p*(1 + s))*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2 + a*p*(1 + s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(4*(-2*c5*s + c4*p*(1 + s))*PhiPrimePrime + a*c4*s*PhiPrimePrimePrime)) - Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*adotoaPrime + 8*(-1)**p*(0,1)**(p/s)*a**2*c4*p**2*(1 + s)**2*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(2*(1 + p + p/s))*adotoaPrime**2 + 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(a*c4*p*(1 + s)*(s + 4*p*(1 + s))*PhiPrimePrime*adotoaPrime + PhiPrime*(2*(s + 2*p*(1 + s))*(-2*c5*s + c4*p*(1 + s))*adotoaPrime + a*c4*p*s*(1 + s)*adotoaPrimePrime)))/((-1)**(p*(2 + 1/s))*(0,1)**(p/s)*s**2*XDS**(p*(2 + 3/(2.*s)))*PhiPrime**2)

          self%EFTc%yp(ind)      = (-4*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*p**2*s**3*XDS**(p + (3*p)/(2.*s))*adotoa**(1 + 2*p)*PhiPrime**(2 + 2*p)*PhiPrimePrime + Sqrt(2)*a*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(2 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*(-3*s*PhiPrime**2 + a**2*(s - p*(1 + 2*s))*PhiPrimePrime**2 + a*PhiPrime*(3*(p + 2*p*s)*PhiPrimePrime - a*s*PhiPrimePrimePrime)) + 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(3 + (2*p*(1 + s))/s)*(2*s*(-28*c5*s**2 + c4*p*(1 + s)*(-s + 6*p*(1 + s)))*PhiPrime**(3 + (2*p*(1 + s))/s) + 2*a**3*c4*p*(1 + s)*(s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**3 + a**2*p*(1 + s)*(-s + 2*p*(1 + s))*PhiPrime**(1 + (2*p*(1 + s))/s)*PhiPrimePrime*((-3*c4*s - 8*c5*s + 4*c4*p*(1 + s))*PhiPrimePrime + 3*a*c4*s*PhiPrimePrimePrime) + a*p*(1 + s)*PhiPrime**(2*(1 + p + p/s))*(-((-64*c5*s**2 + c4*(-3*s**2 + 2*p*s*(1 + s) + 12*p**2*(1 + s)**2))*PhiPrimePrime) + a*s*((-3*c4*s - 8*c5*s + 4*c4*p*(1 + s))*PhiPrimePrimePrime + a*c4*s*PhiPrimePrimePrimePrime))) - 4*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*p**2*s**3*XDS**(p + (3*p)/(2.*s))*adotoa**(2*p)*PhiPrime**(3 + 2*p)*adotoaPrime - Sqrt(2)*a**3*c3*E**((I*p*Pi*(1 + s))/s)*Hds*p*s*(1 + 2*s)*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(3 + p*(2 + 1/s))*adotoaPrime**2 + 16*(-1)**p*(0,1)**(p/s)*a**3*c4*p**3*(1 + s)**3*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(3 + (2*p*(1 + s))/s)*adotoaPrime**3 - Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*(a*(s + p*(2 + 4*s))*PhiPrimePrime*adotoaPrime + PhiPrime*(-3*(p + s + 2*p*s)*adotoaPrime + a*s*adotoaPrimePrime)) + 4*(-1)**p*(0,1)**(p/s)*a**2*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(2*(1 + p + p/s))*adotoaPrime*(a*c4*p*(1 + s)*(s**2 + 6*p*s*(1 + s) + 12*p**2*(1 + s)**2)*PhiPrimePrime*adotoaPrime - c4*p*s**2*PhiPrime*(adotoaPrime - a*adotoaPrimePrime) + PhiPrime*((-4*c5*s*(s + 2*p*(1 + s))**2 + c4*p*(-s**3 + 2*p*s*(1 + s)**2 + 8*p**2*(1 + s)**3))*adotoaPrime + a*c4*p*s*(s**2 + 6*p*(1 + s)**2)*adotoaPrimePrime)) + 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*PhiPrime*(3*a**2*c4*p*(1 + s)*(-s**2 + 4*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2*adotoaPrime +a*p*(1 + s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(3*a*c4*s*(s + 2*p*(1 + s))*PhiPrimePrimePrime*adotoaPrime + PhiPrimePrime*(2*(-4*c5*s*(3*s + 4*p*(1 + s)) + c4*(-3*s**2 + 8*p**2*(1 + s)**2))*adotoaPrime + a*c4*s*(s + 6*p*(1 + s))*adotoaPrimePrime)) - PhiPrime**(2*(1 + p + p/s))*((-4*c5*s**2*(15*s + 16*p*(1 + s)) + c4*p*(1 + s)*(-3*s**2 + 14*p*s*(1 + s) + 12*p**2*(1 + s)**2))*adotoaPrime + a*s*((4*c5*s*(s + 2*p*(1 + s)) - c4*p*(1 + s)*(-s + 4*p*(1 + s)))*adotoaPrimePrime - a*c4*p*s*(1 + s)*adotoaPrimePrimePrime))))/(2.*(-1)**(p*(2 + 1/s))*(0,1)**(p/s)*s**3*XDS**(p*(2 + 3/(2.*s)))*PhiPrime**3)

          self%EFTLambda%yp(ind) =         (-2*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*p*s**3*XDS**(p + (3*p)/(2.*s))*adotoa**(1 + 2*p)*PhiPrime**(2 + 2*p)*PhiPrimePrime - Sqrt(2)*a**3*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(2 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*((p - s + 2*p*s)*PhiPrimePrime**2 + s*PhiPrime*PhiPrimePrimePrime) + 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(3 + (2*p*(1 + s))/s)*(-2*s**2*(-2*c5*s + c4*p*(1 + s))*PhiPrime**(3 + (2*p*(1 + s))/s) + 2*a**3*c4*p*(1 + s)*(s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**3 + a**2*p*(1 + s)*(-s + 2*p*(1 + s))*PhiPrime**(1 + (2*p*(1 + s))/s)*PhiPrimePrime*(4*(-2*c5*s + c4*p*(1 + s))*PhiPrimePrime + 3*a*c4*s*PhiPrimePrimePrime) + a*p*s*(1 + s)*PhiPrime**(2*(1 + p + p/s))*(4*c5*s*PhiPrimePrime - 2*c4*p*(1 + s)*PhiPrimePrime + 4*a*(-2*c5*s + c4*p*(1 + s))*PhiPrimePrimePrime + a**2*c4*s*PhiPrimePrimePrimePrime)) - 2*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*p*s**3*XDS**(p + (3*p)/(2.*s))*adotoa**(2*p)*PhiPrime**(3 + 2*p)*adotoaPrime - Sqrt(2)*a**3*c3*E**((I*p*Pi*(1 + s))/s)*Hds*p*s*(1 + 2*s)*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(3 + p*(2 + 1/s))*adotoaPrime**2 + 16*(-1)**p*(0,1)**(p/s)*a**3*c4*p**3*(1 + s)**3*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(3 + (2*p*(1 + s))/s)*adotoaPrime**3 - Sqrt(2)*a**3*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*((s + p*(2 + 4*s))*PhiPrimePrime*adotoaPrime + s*PhiPrime*adotoaPrimePrime) + 4*(-1)**p*(0,1)**(p/s)*a**2*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**2*adotoaPrime*(a*c4*p*(1 + s)*(s**2 + 6*p*s*(1 + s) + 12*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime*adotoaPrime + a*c4*p*s**2*PhiPrime**(1 + 2*p + (2*p)/s)*adotoaPrimePrime + PhiPrime**(1 + (2*p*(1 + s))/s)*(2*(s + 2*p*(1 + s))**2*(-2*c5*s + c4*p*(1 + s))*adotoaPrime + a*c4*p*s*(s**2 + 6*p*(1 + s)**2)*adotoaPrimePrime)) + 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*PhiPrime*(3*a**2*c4*p*(1 + s)*(-s**2 + 4*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2*adotoaPrime + a*p*(1 + s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(3*a*c4*s*(s + 2*p*(1 + s))*PhiPrimePrimePrime*adotoaPrime + PhiPrimePrime*(4*(3*s + 4*p*(1 + s))*(-2*c5*s + c4*p*(1 + s))*adotoaPrime + a*c4*s*(s + 6*p*(1 + s))*adotoaPrimePrime)) - s*PhiPrime**(2*(1 + p + p/s))*(2*p*(1 + s)*(-2*c5*s + c4*p*(1 + s))*adotoaPrime - a*(2*(s + 2*p*(1 + s))*(-2*c5*s + c4*p*(1 + s))*adotoaPrimePrime + a*c4*p*s*(1 + s)*adotoaPrimePrimePrime))))/((-1)**(p*(2 + 1/s))*(0,1)**(p/s)*s**3*XDS**(p*(2 + 3/(2.*s)))*PhiPrime**3)

          self%EFTgamma1%y(ind)  = ((-4*a**2*c2*Hds**2*(-1 + p)*p*s**3*adotoa**(2*p)*PhiPrime**(2*p))/((-1)**p*XDS**p) + (Sqrt(2)*a*c3*Hds*s*(p - s + 2*p*s)*adotoa**(p*(2 + 1/s))*PhiPrime**(-1 + p*(2 + 1/s))*(adotoa*(3*(p - 2*s + 2*p*s)*PhiPrime + a*s*PhiPrimePrime) + a*s*PhiPrime*adotoaPrime))/((-1)**p*(0,1)**(p/s)*XDS**((p + 2*p*s)/(2.*s))) + (16*c5*s**2*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(-1 + (2*p*(1 + s))/s)*(2*a*(p - 2*s + p*s)*PhiPrimePrime + s*adotoa*(11*PhiPrime + 4*a*PhiPrimePrime) + a*(s + 2*p*(1 + s))*PhiPrime*adotoaPrime))/((-1)**(p*(1 + 1/s))*XDS**(p*(1 + 1/s))) - (4*c4*p*(1 + s)*adotoa**((2*p*(1 + s))/s)*PhiPrime**(-2 + (2*p*(1 + s))/s)*(adotoa**2*(2*(2*s**2 - 9*p*s*(1 + s) + 6*p**2*(1 + s)**2)*PhiPrime**2 + a**2*s*(-s + 2*p*(1 + s))*PhiPrimePrime**2 + a*s*PhiPrime*((-3*s + 4*p*(1 + s))*PhiPrimePrime + a*s*PhiPrimePrimePrime)) + 2*a**2*p*s*(1 + s)*PhiPrime**2*adotoaPrime**2 + a*s*adotoa*PhiPrime*(a*(s + 4*p*(1 + s))*PhiPrimePrime*adotoaPrime + PhiPrime*((-s + 4*p*(1 + s))*adotoaPrime + a*s*adotoaPrimePrime))))/((-1)**(p*(1 + 1/s))*XDS**(p*(1 + 1/s))))/(4.*a**2*H0**2*s**3)

          self%EFTgamma1%yp(ind)  = (-8*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*(-1 + p)*p**2*s**4*XDS**(p + (3*p)/(2.*s))*adotoa**(1 + 2*p)*PhiPrime**(2 + 2*p)*PhiPrimePrime + Sqrt(2)*a*c3*Hds*s*XDS**(p*(1 + 1/s))*adotoa**(2 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*(-3*(-1)**(p*(1 + 1/s))*s*(2*s**2 - 3*p*s*(1 + 2*s) + (p + 2*p*s)**2)*PhiPrime**2 + (-1)**(p*(1 + 1/s))*a**2*s*(p - s + 2*p*s)**2*PhiPrimePrime**2 + a*E**((I*p*Pi*(1 + s))/s)*(p - s + 2*p*s)*PhiPrime*(3*p*(1 + 2*s)*(p - 2*s + 2*p*s)*PhiPrimePrime + a*s**2*PhiPrimePrimePrime)) + 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(3 + (2*p*(1 + s))/s)*(4*s*(-22*c5*s**3 + c4*p*(1 + s)*(2*s**2 - 9*p*s*(1 + s) + 6*p**2*(1 + s)**2))*PhiPrime**(3 + (2*p*(1 + s))/s) - 2*a**3*c4*p*s*(1 + s)*(s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2)*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**3 - a**2*s*(-s + 2*p*(1 + s))*PhiPrime**(1 + (2*p*(1 + s))/s)*PhiPrimePrime*((-16*c5*s**2 + c4*p*(1 + s)*(-3*s + 4*p*(1 + s)))*PhiPrimePrime + 3*a*c4*p*s*(1 + s)*PhiPrimePrimePrime) - a*PhiPrime**(2*(1 + p + p/s))*((-8*c5*s**3*(-2*s + 11*p*(1 + s)) + c4*p*(1 + s)*(3*s**3 + 4*p*s**2*(1 + s) - 36*p**2*s*(1 + s)**2 + 24*p**3*(1 + s)**3))*PhiPrimePrime + a*s**2*((-16*c5*s**2 + c4*p*(1 + s)*(-3*s + 4*p*(1 + s)))*PhiPrimePrimePrime + a*c4*p*s*(1 + s)*PhiPrimePrimePrimePrime))) - 8*a**3*c2*E**((I*p*Pi*(3 + 2*s))/(2.*s))*Hds**2*(-1 + p)*p**2*s**4*XDS**(p + (3*p)/(2.*s))*adotoa**(2*p)*PhiPrime**(3 + 2*p)*adotoaPrime + Sqrt(2)*a**3*c3*E**((I*p*Pi*(1 + s))/s)*Hds*p*s**2*(1 + 2*s)*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(3 + p*(2 + 1/s))*adotoaPrime**2 - 16*(-1)**p*(0,1)**(p/s)*a**3*c4*p**3*s*(1 + s)**3*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(3 + (2*p*(1 + s))/s)*adotoaPrime**3 + Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*(a*s*(s + p*(2 + 4*s))*PhiPrimePrime*adotoaPrime + PhiPrime*(3*(-2*s**2 - p*s*(1 + 2*s) + (p + 2*p*s)**2)*adotoaPrime + a*s**2*adotoaPrimePrime)) - 4*(-1)**p*(0,1)**(p/s)*a**2*s*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(2*(1 + p + p/s))*adotoaPrime*(PhiPrimePrime*(-8*c5*s*(-2*s**2 - 3*p*s*(1 + s) + 2*p**2*(1 + s)**2) + a*c4*p*(1 + s)*(s**2 + 6*p*s*(1 + s) + 12*p**2*(1 + s)**2)*adotoaPrime) - c4*p*s**2*PhiPrime*(adotoaPrime - a*adotoaPrimePrime) + PhiPrime*((-4*c5*s*(s + 2*p*(1 + s))**2 + c4*p*(-s**3 + 2*p*s*(1 + s)**2 + 8*p**2*(1 + s)**3))*adotoaPrime + a*c4*p*s*(s**2 + 6*p*(1 + s)**2)*adotoaPrimePrime)) - 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*PhiPrime*(a*s*(-s + 2*p*(1 + s))*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2*(-8*c5*s*(p - 2*s + p*s) + 3*a*c4*p*(1 + s)*(s + 2*p*(1 + s))*adotoaPrime) + s*PhiPrime**(1 + (2*p*(1 + s))/s)*(a*s*PhiPrimePrimePrime*(-8*c5*s*(p - 2*s + p*s) + 3*a*c4*p*(1 + s)*(s + 2*p*(1 + s))*adotoaPrime) + PhiPrimePrime*(2*a*(-4*c5*s*(4*s**2 + 5*p*s*(1 + s) + 2*p**2*(1 + s)**2) + c4*p*(1 + s)*(-3*s**2 + 8*p**2*(1 + s)**2))*adotoaPrime + s*(8*c5*s*(p - 2*s + p*s) + a**2*c4*p*(1 + s)*(s + 6*p*(1 + s))*adotoaPrimePrime))) + PhiPrime**(2*(1 + p + p/s))*((-4*c5*s**3*(21*s + 20*p*(1 + s)) + c4*p*(1 + s)*(9*s**3 - 32*p*s**2*(1 + s) - 12*p**2*s*(1 + s)**2 + 24*p**3*(1 + s)**3))*adotoaPrime + a*s**2*((-4*c5*s*(s + 2*p*(1 + s)) + c4*p*(1 + s)*(-s + 4*p*(1 + s)))*adotoaPrimePrime + a*c4*p*s*(1 + s)*adotoaPrimePrimePrime))))/(4.*(-1)**(p*(2 + 1/s))*(0,1)**(p/s)*a**3*H0**2*s**4*XDS**(p*(2 + 3/(2.*s)))*adotoa*PhiPrime**3)


          self%EFTgamma2%y(ind)  =  (-(Sqrt(2)*a*c3*E**((I*p*Pi*(1 + s))/s)*Hds*s*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))) + 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**((2*p*(1 + s))/s)*(2*(-8*c5*s**2 + c4*p*(1 + s)*(-s + 2*p*(1 + s)))*PhiPrime + a*c4*p*s*(1 + s)*PhiPrimePrime) + 4*(-1)**p*(0,1)**(p/s)*a*c4*p*s*(1 + s)*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(1 + (2*p*(1 + s))/s)*adotoaPrime)/((-1)**(p*(2 + 1/s))*(0,1)**(p/s)*a*H0*s**2*XDS**(p*(2 + 3/(2.*s)))*PhiPrime)

          self%EFTgamma2%yp(ind)  =  (-(Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*p*s*(1 + 2*s)*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(1 + p*(2 + 1/s))*PhiPrime**(1 + p*(2 + 1/s))*PhiPrimePrime) - 4*(-1)**p*(0,1)**(p/s)*XDS**(p + p/(2.*s))*adotoa**(2*(1 + p + p/s))*(2*s*(-8*c5*s**2 + c4*p*(1 + s)*(-s + 2*p*(1 + s)))*PhiPrime**(2*(1 + p + p/s)) - a**2*c4*p*s*(1 + s)*(-s + 2*p*(1 + s))*PhiPrime**((2*p*(1 + s))/s)*PhiPrimePrime**2 - a*p*(1 + s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(4*(-8*c5*s**2 + c4*p*(1 + s)*(-s + 2*p*(1 + s)))*PhiPrimePrime + a*c4*s**2*PhiPrimePrimePrime)) - Sqrt(2)*a**2*c3*E**((I*p*Pi*(1 + s))/s)*Hds*p*s*(1 + 2*s)*(p - s + 2*p*s)*XDS**(p*(1 + 1/s))*adotoa**(p*(2 + 1/s))*PhiPrime**(2 + p*(2 + 1/s))*adotoaPrime + 8*(-1)**p*(0,1)**(p/s)*a**2*c4*p**2*s*(1 + s)**2*XDS**(p + p/(2.*s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**(2*(1 + p + p/s))*adotoaPrime**2 + 4*(-1)**p*(0,1)**(p/s)*a*XDS**(p + p/(2.*s))*adotoa**(1 + (2*p*(1 + s))/s)*PhiPrime**(1 + (2*p*(1 + s))/s)*(a*c4*p*s*(1 + s)*(s + 4*p*(1 + s))*PhiPrimePrime*adotoaPrime + PhiPrime*(2*(s + 2*p*(1 + s))*(-8*c5*s**2 + c4*p*(1 + s)*(-s + 2*p*(1 + s)))*adotoaPrime + a*c4*p*s**2*(1 + s)*adotoaPrimePrime)))/((-1)**(p*(2 + 1/s))*(0,1)**(p/s)*a**2*H0*s**3*XDS**(p*(2 + 3/(2.*s)))*adotoa*PhiPrime**2)

          self%EFTgamma3%y(ind)  = (4*(-2*c5*s + c4*p*(1 + s))*adotoa**((2*p*(1 + s))/s)*PhiPrime**((2*p*(1 + s))/s))/((-1)**(p*(1 + 1/s))*s*XDS**(p*(1 + 1/s)))

          self%EFTgamma3%yp(ind)  =  (8*p*(1 + s)*(-2*c5*s + c4*p*(1 + s))*adotoa**(-1 + (2*p*(1 + s))/s)*PhiPrime**(-1 + (2*p*(1 + s))/s)*(adotoa*PhiPrimePrime + PhiPrime*adotoaPrime))/((-1)**(p*(1 + 1/s))*s**2*XDS**(p*(1 + 1/s)))

          self%EFTgamma4%y(ind)  = -self%EFTgamma3%y(ind)

          self%EFTgamma4%yp(ind)  = -self%EFTgamma3%yp(ind)

          self%EFTgamma4%ypp(ind) = (-8*p*(1 + s)*(-2*c5*s + c4*p*(1 + s))*adotoa**(-2 + (2*p*(1 + s))/s)*PhiPrime**(-2 + (2*p*(1 + s))/s)*(adotoa**2*((-s + 2*p*(1 + s))*PhiPrimePrime**2 + s*PhiPrime*PhiPrimePrimePrime) + (-s + 2*p*(1 + s))*PhiPrime**2*adotoaPrime**2 + adotoa*PhiPrime*(4*p*(1 + s)*PhiPrimePrime*adotoaPrime + s*PhiPrime*adotoaPrimePrime)))/((-1)**(p*(1 + 1/s))*s**3*XDS**(p*(1 + 1/s)))

         ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
         ! modification ends
         ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          end subroutine

    end subroutine EFTCAMBWBGalileonSolveBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBWBGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_WB_Galileon)  :: self   !< the base class

        self%parameter_number = 3 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    end subroutine EFTCAMBWBGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBWBGalileonFeedback( self, print_params )

        implicit none

        class(EFTCAMB_WB_Galileon)    :: self         !< the base class
        logical, optional                  :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        logical                            :: print_params_temp

        ! print general model informations:
        if ( self%WBGalileon_c2 == -1._dl ) then

          write(*,*)
          write(*,'(a,a)')    '   Model               =  ', self%name
          write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
          write(*,'(a,F12.6)')   '                 WB_Galileon_XDS    ='  , self%XDS
          write(*,'(a,F12.6)')   '                 WB_Galileon_p    ='  , self%p
          write(*,'(a,F12.6)')   '                 WB_Galileon_q    ='  , self%q
          write(*,'(a,F12.6)')   '                 c2    ='  , self%c2
          write(*,'(a,F12.6)')   '                 c3    ='  , self%c3
          write(*,'(a,F12.6)')   '                 c4    ='  , self%c4
          write(*,'(a,F12.6)')   '                 c5    ='  , self%c5

        end if
        ! print the values of the parameters:
        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

    end subroutine EFTCAMBWBGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBWBGalileonParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_WB_Galileon)   :: self   !< the base class
        integer     , intent(in)          :: i      !< the index of the parameter
        character(*), intent(out)         :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            name = 'WB_Galileon_XDS'   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            return
        end if
        if ( i==2 ) then               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            name = 'WB_Galileon_p'     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            return                     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        end if                         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if ( i==3 ) then               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            name = 'WB_Galileon_s'     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            return                     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        end if                         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBWBGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBWBGalileonParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_WB_Galileon)   :: self       !< the base class
        integer     , intent(in)          :: i          !< The index of the parameter
        character(*), intent(out)         :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==1 ) then
            latexname = '\xi'
            return
        end if
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBWBGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBWBGalileonParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_WB_Galileon)   :: self   !< the base class
        integer , intent(in)              :: i      !< The index of the parameter
        real(dl), intent(out)             :: value  !< the output value of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_value.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==0 ) then
            value = 0._dl
            return
        end if
        if ( i==1 ) then
            value = self%csi
            return
        end if

    end subroutine EFTCAMBWBGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBWBGalileonBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        if(x>=self%x_final) return
        if(x<=self%x_initial) return

        if(a==0._dl) then
            return
        else if (eft_cache%adotoa==0._dl) then
            call self%compute_adotoa( a, eft_par_cache, eft_cache )
            call self%compute_H_derivs( a, eft_par_cache, eft_cache )
            if(eft_cache%adotoa==0._dl) return
            if(eft_cache%Hdot==0._dl) return
            if(eft_cache%Hdotdot==0._dl) return
        end if

        call self%EFTOmega%precompute(x, ind, mu )

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = self%EFTOmega%value( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaP    = self%EFTOmega%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPP   = self%EFTOmega%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPPP  = self%EFTOmega%third_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBWBGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBWBGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        if(x>=self%x_final) return
        if(x<=self%x_initial) return

        call self%EFTgamma1%precompute(x, ind, mu )
        !
        ! ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = self%EFTgamma1%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma1P  = self%EFTgamma1%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma2V  = self%EFTgamma2%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma2P  = self%EFTgamma2%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma3V  = self%EFTgamma3%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma3P  = self%EFTgamma3%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4V  = self%EFTgamma4%value( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4P  = self%EFTgamma4%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma4PP = self%EFTgamma4%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTGamma5V  = 0.5_dl*eft_cache%EFTGamma3V
        eft_cache%EFTGamma5P  = 0.5_dl*eft_cache%EFTGamma3P
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBWBGalileonSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBWBGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot

        a2 = a*a
        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
        temp = 0.5_dl*a2*(eft_par_cache%h0_Mpc)**2*( Omega_tot + sqrt( Omega_tot**2 +4._dl*eft_par_cache%omegav ) )
        eft_cache%adotoa = sqrt( temp )

    end subroutine EFTCAMBWBGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBWBGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime,Omega_tot_primeprimeprime, Omega_phi0,Omega_tot_primeprimeprimeprime
        !
        a2 = a*a
        !
        if(a*eft_cache%adotoa==0._dl) return

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
        Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                          & -(eft_cache%grhonu_tot+eft_cache%gpinu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)
        Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                          & +(4._dl*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)-eft_cache%gpinudot_tot/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)
        Omega_phi0 = eft_par_cache%omegav
        eft_cache%Hdot = eft_cache%adotoa**2 +0.25_dl*(eft_par_cache%h0_Mpc)**2*a**3*( 1._dl + Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_prime
        eft_cache%Hdotdot = 2._dl*eft_cache%adotoa*eft_cache%Hdot +3._dl*eft_cache%adotoa*( eft_cache%Hdot -eft_cache%adotoa**2 ) +0.25_dl*(eft_par_cache%h0_Mpc)**2*eft_cache%adotoa*a2**2&
            & *( ( 1._dl +Omega_tot/sqrt( Omega_tot**2 +4._dl*Omega_phi0 ) )*Omega_tot_primeprime +Omega_tot_prime**2&
            & *( 4._dl*Omega_phi0/( Omega_tot**2 +4._dl*Omega_phi0 )**( 1.5_dl ) ) )

    end subroutine EFTCAMBWBGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBWBGalileonAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_WB_Galileon)              :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBWBGalileonAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBWBGalileonAdditionalModelStability = .True.

    end function EFTCAMBWBGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_WB_Galileon
