module nonlinear_transfer
use vars
use SEE
use rt_coef
use maths
use ogpf
implicit none
contains

!------------------------------------------------------------
! Do a synthesis calling the appropriate routines
!------------------------------------------------------------
    subroutine do_transfer( in_params, in_fixed, in_observation, slab, output, error )
    type(variable_parameters) :: in_params
    type(type_observation) :: in_observation
    type(fixed_parameters) :: in_fixed
    type(type_slab) :: slab
    real(kind=8) :: output(0:3,in_fixed%no)
    integer :: i, iteration, loop_shell, error, j
    real(kind=8) :: I0, Q0, U0, V0, ds, Imax, Ic, factor, eta0
    real(kind=8) :: tolerance(2), vmacro
    real(kind=8), allocatable ::          &
      epsI(:), epsQ(:), epsU(:), epsV(:), &
      etaI(:), etaQ(:), etaU(:), etaV(:), &
               rhoQ(:), rhoU(:), rhoV(:)  !
    real(kind=8), allocatable, dimension(:) :: dtau, delta, prof(:)
    real(kind=8), allocatable :: m1(:,:), m2(:,:), Stokes0(:)
    real(kind=8), allocatable :: O_evol(:,:), psi_matrix(:,:), J00(:), J20(:), J00_nu(:,:), J20_nu(:,:)

    integer(kind=4) :: &
      line, layer, angle, spectrum_size, begin_f, end_f, n_points, f, save_nemiss, k, km, kp, kto, kfrom, kstep
    real(kind=8) :: &
      mu, illumination_cone_cosine, illumination_cone_sine, cos_alpha, &
      wavelength, v_los, d_nu_dopp, resolution, adamp, save_thetad, save_chid, save_gammad, &
      eta_I

    ! A ray is considered horizontal it its mu is less than
    real(kind=8), parameter :: horizontal_mu_limit = 0.01

    logical :: is_horizontal_ray

    ! Quantities for the formal solution.
    real(kind=8), allocatable :: &
      IQUV( :, : ), &
      chim( : ), chi0( : ), chip( : ), dtm( : ), dtp( : ), exu( : ), &
      sm( :, : ), s0( :, : ), sp( :, : ), &
      ab_matrix( :, :, :, : ), source_vector( :, :, : ), &
      tau_nu( : ), tau( : )
    real(kind=8) :: dm, dp, psim, psi0, psip, mat1( 4, 4 ), mat2( 4, 4 )

    ! ogpf
    type(gpf) :: gp
    character(len=4) :: tmp_str

      if ( verbose_mode == 1 ) then
          print *, 'Starting transfer...'
      endif

      spectrum_size = multiplets( atom%ntran )%end

      ! Generate the frequency and wavelength grids for all multiplets.
      do line = 1, atom%ntran
        n_points   = multiplets( line )%no
        resolution = ( multiplets( line )%d_wl_max - multiplets( line )%d_wl_min ) / dble( n_points - 1 )
        wavelength = multiplets( line )%wl

        allocate( multiplets( line )%lambdas(   n_points ), source = 0d0 )
        allocate( multiplets( line )%nus(       n_points ), source = 0d0 )
        allocate( multiplets( line )%d_lambdas( n_points ), source = 0d0 )
        allocate( multiplets( line )%d_nus(     n_points ), source = 0d0 )

        ! Set an equidistant wavelength grid.
        multiplets( line )%d_lambdas( 1:n_points ) = [ (f, f = 0, n_points - 1) ] * resolution + multiplets( line )%d_wl_min
        multiplets( line )%lambdas(   1:n_points ) = multiplets( line )%d_lambdas( 1:n_points ) + wavelength

        ! Convert Ångströms to cm with a factor of 10^{-8} to get Hz.
        multiplets( line )%nus( 1:n_points )   = ( PC / 1d-8 ) / multiplets( line )%lambdas( 1:n_points )
        multiplets( line )%d_nus( 1:n_points ) = multiplets( line )%nus( 1:n_points ) - ( PC / 1d-8 / wavelength )
      enddo


      ! Create arrays for the current and the previous radiation field
      ! parameters, \( \bar{n} \) and \( \omega \), ...
      allocate( slab%nbar(      slab%n_layers, atom%ntran ), source = 0d0 )
      allocate( slab%omega(     slab%n_layers, atom%ntran ), source = 0d0 )
      allocate( slab%nbar_old(  slab%n_layers, atom%ntran ), source = 0d0 )
      allocate( slab%omega_old( slab%n_layers, atom%ntran ), source = 0d0 )

      ! ...and initialize the current nbar and omega in all layers with the same
      ! values using Allen's tables to start iterating.
      if ( verbose_mode == 1 ) then
        write (*, '(a)') 'Initial nbar and omega:'
      endif
      do line = 1, atom%ntran
        wavelength = atom%wavelength( line )
        ! Same value for different layers.
        slab%nbar(  :, line ) = nbar_allen(  wavelength, in_fixed, in_params, atom%reduction_factor(       line ) )
        slab%omega( :, line ) = omega_allen( wavelength, in_fixed, in_params, atom%reduction_factor_omega( line ) )
        if ( verbose_mode == 1 ) then
          write (*, '(4x, a, f9.3, 2(a, e8.2))') 'wavelength = ', wavelength, ' A, nbar = ', slab%nbar( 1, line ), ', omega = ', slab%omega( 1, line )
        endif
      enddo
      ! Set the previous values to the current ones.
      slab%nbar_old  = slab%nbar
      slab%omega_old = slab%omega


      ! Create an array for the boundary conditions...
      allocate( slab%boundary( 4, spectrum_size, slab%aq_size ), source = 0d0 )

      ! From Eq. (12.32) of Landi degl'Innocenti & Landolfi (2004) obtain
      ! properties of the illumination cone for the given height...
      illumination_cone_sine   = RSUN / ( RSUN + in_params%height )
      illumination_cone_cosine = sqrt( 1d0 - illumination_cone_sine**2 )
      if ( verbose_mode == 1 ) then
          write (*, '(a, f5.2)') 'Illumination cone mu = ', illumination_cone_cosine
      endif

      ! ...and use Allen's tables for the center-to-limb variation of the
      ! continuum intensity to set all Stokes I values for the current angle.
      if ( verbose_mode == 1 ) then
          write (*, '(a)') 'Boundary intensities:'
      endif
      do angle = 1, slab%aq_size
        ! for all inclinations \( \theta \) having
        ! \( \mu = \cos\theta > \cos\gamma \).
        mu = cos( slab%aq_inclination( angle ) * PI / 180d0 )
        if ( mu > illumination_cone_cosine ) then
          ! Eq. 12.33 of Landi degl'Innocenti & Landolfi (2004).
          cos_alpha = sqrt( mu**2 - illumination_cone_cosine**2 ) / illumination_cone_sine
          do line = 1, atom%ntran
            wavelength = multiplets( line )%wl
            begin_f    = multiplets( line )%begin
            end_f      = multiplets( line )%end
            slab%boundary( 1, begin_f:end_f, angle ) = I0_allen( wavelength, cos_alpha )
          enddo
        else
          ! For other angles there is no incoming boundary intensity and it is
          ! already set to zero.
          continue
        endif
        if ( verbose_mode == 1 ) then
          write (*, '(4x, a, f6.3, a, *(e9.2, x))') 'mu = ', mu, ', I0 = ', slab%boundary( 1, multiplets%begin, angle )
        endif
      enddo


      ! Generate the absorption profile.
      allocate( slab%absorption_profile(      slab%n_layers, spectrum_size, slab%aq_size ), source = 0d0 )
      allocate( slab%absorption_profile_norm( slab%n_layers,    atom%ntran, slab%aq_size ), source = 0d0 )

      do angle = 1, slab%aq_size
        do layer = 1, slab%n_layers
          ! The LoS velocity is the projection of v_z along the angle direction.
          ! Take the opposite sign following the astrophysical convention that
          ! the velocity is positive for redshifted profiles.  For more details
          ! see a comment below.
          v_los = -slab%v_z( layer ) * cos( slab%aq_inclination( angle ) * PI / 180d0 )
          ! The damping parameter.
          adamp = slab%damping( layer )
          do line = 1, atom%ntran
            wavelength = multiplets( line )%wl
            begin_f    = multiplets( line )%begin
            end_f      = multiplets( line )%end

            ! The Doppler width in Hz.
            d_nu_dopp = slab%vthermal( layer ) * 1d5 / ( wavelength * 1d-8 )

            ! Calculate the profile
            slab%absorption_profile( layer, begin_f:end_f, angle ) = dble( &
                profile( adamp, -multiplets( line )%d_nus / d_nu_dopp - v_los / slab%vthermal( layer ) ) &
                / ( d_nu_dopp * SQRTPI ) &
              )

            ! Calculate the norm of the (incomplete) profile.  Negate the
            ! frequency shfits as in the calculation of the profile.
            slab%absorption_profile_norm( layer, line, angle ) = &
              int_tabulated( -multiplets( line )%d_nus, slab%absorption_profile( layer, begin_f:end_f, angle ) )
          end do ! 1 <= line <= atom%ntran
        end do ! 1 <= layer <= slab%n_layers
      end do ! 1 <= angle <= slab%aq_size

      !call gp%title( 'Absorption profile along different rays' )
      !call gp%xlabel ( 'd_lambda, A' )
      !call gp%ylabel ( 'phi, [Hz^{-1}]' )
      !call gp%options( 'set style data linespoints; set grid; set key top left' )
      !call gp%plot( [(i * 1d0, i = 1, 200)], slab%absorption_profile( 1, 1:200, 5:6 ) ) ! 'lt 7'
      !stop

      allocate( slab%propagation_matrix( 4, 4, slab%n_layers, spectrum_size, slab%aq_size ), source = 0d0 )
      allocate( slab%emission_vector(       4, slab%n_layers, spectrum_size, slab%aq_size ), source = 0d0 )

      allocate( ab_matrix(     4, 4, slab%n_layers, spectrum_size ), source = 0d0 )
      allocate( source_vector(    4, slab%n_layers, spectrum_size ), source = 0d0 )
      allocate( IQUV(             4,                spectrum_size ), source = 0d0 )

      allocate( chim(  spectrum_size ), source = 0d0 )
      allocate( chi0(  spectrum_size ), source = 0d0 )
      allocate( chip(  spectrum_size ), source = 0d0 )
      allocate( dtm(   spectrum_size ), source = 0d0 )
      allocate( dtp(   spectrum_size ), source = 0d0 )
      allocate( exu(   spectrum_size ), source = 0d0 )
      allocate( sm( 4, spectrum_size ), source = 0d0 )
      allocate( s0( 4, spectrum_size ), source = 0d0 )
      allocate( sp( 4, spectrum_size ), source = 0d0 )

      allocate( tau_nu( spectrum_size ), source = 0d0 )
      allocate( tau( atom%ntran ), source = 0d0 )

      !allocate(J00_nu(slab%n_layers,in_fixed%no))
      !allocate(J20_nu(slab%n_layers,in_fixed%no))
      !allocate(J00(slab%n_layers))
      !allocate(J20(slab%n_layers))

      ! Initialize iterations
      tolerance = 100.d0
      iteration = 1

      do while ( iteration < 50 .and. maxval( tolerance ) > 1d-3 )
        ! Iterate along the layers computing and storing the radiative trasnfer
        ! quantities because fill_SEE() evaluates the density matrix elements
        ! only for the current layer and is very time-consuming.
        do layer = 1, slab%n_layers

          ! To fill and solve the statistical equilibrium equations, get for
          ! the current slab the following parameters needed by fill_SEE().  The
          ! original values provided in init_parameters.dat are not needed.
          !in_params%delta_collision  ! Don't change this.
          !in_params%bgauss2          ! Don't use the 2nd component in the non-linear mode.
          !in_params%chibd2
          !in_params%thetabd2
          in_params%bgauss  = slab%B(    layer )
          in_params%thetabd = slab%thB(  layer )
          in_params%chibd   = slab%chiB( layer )

          ! ...pass the new nbar and omega in the current layer.
          nbarExternal  = slab%nbar(  layer, : )
          omegaExternal = slab%omega( layer, : )


          ! For the first component only...
          call fill_SEE( in_params, in_fixed, 1, error )

          ! To calculate the absorption/emission coefficients, set the line
          ! broadening parameters...
          in_params%damping = slab%damping(  layer )
          in_params%vdopp   = slab%vthermal( layer )

          ! Save the original line number for which the intensity output must be
          ! computed.
          save_nemiss = in_fixed%nemiss

          ! To save time on allocating/deallocating arrays, first iterate lines.
          do line = 1, atom%ntran

            ! Set the current line properties.
            in_fixed%nemiss = line

            in_fixed%no       = multiplets( line )%no
            begin_f           = multiplets( line )%begin
            end_f             = multiplets( line )%end
            in_fixed%wl       = multiplets( line )%wl
            in_fixed%d_wl_min = multiplets( line )%d_wl_min
            in_fixed%d_wl_max = multiplets( line )%d_wl_max
            in_observation%wl = multiplets( line )%d_lambdas

            ! TODO: replace this with static arrays when possible.
            if ( .not.allocated( epsI ) ) allocate( epsI( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( epsQ ) ) allocate( epsQ( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( epsU ) ) allocate( epsU( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( epsV ) ) allocate( epsV( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( etaI ) ) allocate( etaI( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( etaQ ) ) allocate( etaQ( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( etaU ) ) allocate( etaU( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( etaV ) ) allocate( etaV( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( rhoQ ) ) allocate( rhoQ( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( rhoU ) ) allocate( rhoU( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( rhoV ) ) allocate( rhoV( in_fixed%no ), source = 0d0 )
            if ( .not.allocated( dtau ) ) allocate( dtau( in_fixed%no ), source = 0d0 )

            ! Save the original LoS.
            save_thetad = in_fixed%thetad
            save_chid   = in_fixed%chid
            save_gammad = in_fixed%gammad

            do angle = 1, slab%aq_size

              ! Set the LoS angles to the current ray of the quadrature.
              in_fixed%thetad = slab%aq_inclination( angle )
              in_fixed%chid   = slab%aq_azimuth(     angle )
              !in_fixed%gammad = in_fixed%gammad ! TODO: not clear what should be here.

              ! Set the LoS velocity: the reversed projection of v_z along the
              ! current LoS.
              in_params%vmacro = -slab%v_z( layer ) * cos( in_fixed%thetad * PI / 180d0 )

              call calc_rt_coef( in_params, in_fixed, in_observation, 1 )

              ! We have the following arrays of the ( 0:3, n_nu ) size:
              !  epsilon       epsilon_zeeman
              !  eta           eta_zeeman
              !  eta_stim      eta_stim_zeeman
              !  mag_opt       mag_opt_zeeman
              !  mag_opt_stim  mag_opt_stim_zeeman

              ! Collect the elements of the absorption matrix and the emision vector.
              if ( abs( in_fixed%use_atomic_pol ) == 1 ) then
                ! Emission
                epsI = epsilon( 0, : )
                epsQ = epsilon( 1, : )
                epsU = epsilon( 2, : )
                epsV = epsilon( 3, : )
                ! Absorption including stimulated emission
                etaI = eta( 0, : ) - use_stim_emission_RT * eta_stim( 0, : ) 
                etaQ = eta( 1, : ) - use_stim_emission_RT * eta_stim( 1, : )
                etaU = eta( 2, : ) - use_stim_emission_RT * eta_stim( 2, : )
                etaV = eta( 3, : ) - use_stim_emission_RT * eta_stim( 3, : )
                ! Magneto-optical effects
                if ( use_mag_opt_RT == 1 ) then
                  rhoQ = mag_opt( 1, : ) - use_stim_emission_RT * mag_opt_stim( 1, : )
                  rhoU = mag_opt( 2, : ) - use_stim_emission_RT * mag_opt_stim( 2, : )
                  rhoV = mag_opt( 3, : ) - use_stim_emission_RT * mag_opt_stim( 3, : )
                else
                  rhoQ = 0d0
                  rhoU = 0d0
                  rhoV = 0d0
                endif
              else ! No atomic polarization, the Zeeman effect only.
                ! Emission
                epsI = epsilon_zeeman( 0, : )
                epsQ = epsilon_zeeman( 1, : )
                epsU = epsilon_zeeman( 2, : )
                epsV = epsilon_zeeman( 3, : )
                ! Absorption including stimulated emission
                etaI = eta_zeeman( 0, : ) - use_stim_emission_RT * eta_stim_zeeman( 0, : ) + 1d-20
                etaQ = eta_zeeman( 1, : ) - use_stim_emission_RT * eta_stim_zeeman( 1, : ) + 1d-20
                etaU = eta_zeeman( 2, : ) - use_stim_emission_RT * eta_stim_zeeman( 2, : ) + 1d-20
                etaV = eta_zeeman( 3, : ) - use_stim_emission_RT * eta_stim_zeeman( 3, : ) + 1d-20
                ! Magneto-optical terms
                if ( use_mag_opt_RT == 1 ) then
                  rhoQ = mag_opt_zeeman( 1, : ) - use_stim_emission_RT * mag_opt_stim_zeeman( 1, : )
                  rhoU = mag_opt_zeeman( 2, : ) - use_stim_emission_RT * mag_opt_stim_zeeman( 2, : )
                  rhoV = mag_opt_zeeman( 3, : ) - use_stim_emission_RT * mag_opt_stim_zeeman( 3, : )
                else
                  rhoQ = 0d0
                  rhoU = 0d0
                  rhoV = 0d0
                endif
              endif ! use_atomic_pol

              ! TODO: Replace this with the inline code, subroutine calls are slow for each frequency.
              do f = 1, in_fixed%no
                call fill_absorption_matrix( slab%propagation_matrix( :, :, layer, begin_f + f - 1, angle ), etaI( f ), etaQ( f ), etaU( f ), etaV( f ), rhoQ( f ), rhoU( f ), rhoV( f ) )
              enddo
              slab%emission_vector( 1, layer, begin_f:end_f, angle ) = epsI
              slab%emission_vector( 2, layer, begin_f:end_f, angle ) = epsQ
              slab%emission_vector( 3, layer, begin_f:end_f, angle ) = epsU
              slab%emission_vector( 4, layer, begin_f:end_f, angle ) = epsV

            enddo ! angle

            deallocate( epsI, epsQ, epsU, epsV, etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV, dtau )

            ! Restore the original LoS.
            in_fixed%thetad = save_thetad
            in_fixed%chid   = save_chid
            in_fixed%gammad = save_gammad

          enddo ! line

          ! Scale the absorption matrix and the emission vector to physical
          ! units multiplying them by the number density.
          slab%propagation_matrix( :, :, layer, :, : ) = slab%propagation_matrix( :, :, layer, :, : ) * slab%density( layer )
          slab%emission_vector(       :, layer, :, : ) = slab%emission_vector(       :, layer, :, : ) * slab%density( layer )

          ! Restore the original properties of the line for which the intensity
          ! output must be computed.
          in_fixed%nemiss   = save_nemiss
          in_fixed%no       = multiplets( save_nemiss )%no
          in_fixed%wl       = multiplets( save_nemiss )%wl
          in_fixed%d_wl_min = multiplets( save_nemiss )%d_wl_min
          in_fixed%d_wl_max = multiplets( save_nemiss )%d_wl_max

        enddo ! layer

        ! Now solve the transfer equation and intergrate the radiation field tensors.
        !=======================================================================
        ! First, iterate angles as the last one will be the required LoS for the
        ! intensity output.
        do angle = 1, slab%aq_size

          mu = cos( slab%aq_inclination( angle ) * PI / 180d0 )

          is_horizontal_ray = .false.
          if ( abs( mu ) <= horizontal_mu_limit ) then
            is_horizontal_ray = .true.
            ! The order doesn't matter here as all layers are not connected.
            kfrom = 2
            kto   = slab%n_layers
            kstep = 1
          else if ( mu > horizontal_mu_limit ) then
            ! Going upwards.
            kfrom = 2
            kto   = slab%n_layers
            kstep = 1
          else ! mu < horizontal_mu_limit
            ! Going downwards.
            kfrom = slab%n_layers - 1
            kto   = 1
            kstep = -1
          end if

          ! Initialize with the boundary conditions.
          IQUV( :, : ) = slab%boundary( :, :, angle )

          ! Nullify the optical depth arrays.
          tau_nu( : ) = 0d0
          tau( : )    = 0d0

          ! Calculate J00 and J20
          continue

          ! Slice the absorption matrix and the emission vector for one angle
          ! and reduce them both.
          ab_matrix(     :, :, :, : ) = slab%propagation_matrix( :, :, :, :, angle )
          source_vector(    :, :, : ) = slab%emission_vector(       :, :, :, angle )
          do f = 1, spectrum_size
            do k = 1, slab%n_layers
              eta_I = ab_matrix( 1, 1, k, f )
              ! Normalize the absorption matrix \( \hat{K} \) by \( \eta_I \) to
              ! get \( \hat{K}^* \) and remove the principal diagonal to get
              ! \( \hat{K}^\prime \).
              ab_matrix(  :, :, k, f ) = ab_matrix(  :, :, k, f ) / eta_I - identity_4x4
              ! Normalize the emission vector by \( \eta_I \) as well.
              source_vector( :, k, f ) = source_vector( :, k, f ) / eta_I
            end do ! 1 <= k <= slab%n_layers
          end do ! 1 <= f <= spectrum_size

          ! Integrate the transfer equation along all slabs using the method of
          ! short characteristics of the DELO/DELOPAR kind with linear or
          ! parabolic interpolation of the source function.
          do k = kfrom, kto, kstep

            if ( is_horizontal_ray ) then
              ! Asymptotic solution for the horizontal (infinite) ray:
              ! \( \vec{I} = \hat{K}^{-1} \vec{\epsilon} \).
              do f = 1, spectrum_size
                mat1 = slab%propagation_matrix( :, :, k, f, angle )
                call invert( mat1 )
                IQUV( :, f ) = matmul( mat1, slab%emission_vector( :, k, f, angle ) )
              end do
            else
              ! Full formal solution.
              if ( k == kto ) then
                ! Linear short-characteristics at the boundary.
                km = k - kstep
                chim(  : ) = slab%propagation_matrix( 1, 1, km, :, angle )
                chi0(  : ) = slab%propagation_matrix( 1, 1, k,  :, angle )
                chip(  : ) = 0d0
                sm( :, : ) = source_vector( :, km, : )
                s0( :, : ) = source_vector( :, k,  : )
                sp( :, : ) = 0d0
                dm = abs( ( slab%z( k ) - slab%z( km ) ) / mu )
                dp = 0d0
              else
                ! Parabolic short-characteristics inbetween.
                km = k - kstep
                kp = k + kstep
                chim(  : ) = slab%propagation_matrix( 1, 1, km, :, angle )
                chi0(  : ) = slab%propagation_matrix( 1, 1, k,  :, angle )
                chip(  : ) = slab%propagation_matrix( 1, 1, kp, :, angle )
                sm( :, : ) = source_vector( :, km, : )
                s0( :, : ) = source_vector( :, k,  : )
                sp( :, : ) = source_vector( :, kp, : )
                dm = abs( ( slab%z( k  ) - slab%z( km ) ) / mu )
                dp = abs( ( slab%z( kp ) - slab%z( k  ) ) / mu )
              end if

              dtm = 0.5d0 * dm * ( chim + chi0 ) 
              dtp = 0.5d0 * dp * ( chi0 + chip ) 

              ! Increment the optical depth with the previous (km->k0) optical depth element.
              tau_nu = tau_nu + dtm

              where ( dtm >= 1d-4 )
                exu = exp( -dtm )
              else where
                exu = ( dtm * 0.5d0 - 1d0 ) * dtm + 1d0
              end where

              do f = 1, spectrum_size

                call lin_sc( dtm( f ), psim, psi0 )
                mat1 = exu( f ) * identity_4x4 - psim * ab_matrix( :, :, km, f )
                mat2 =            identity_4x4 + psi0 * ab_matrix( :, :, k,  f )
                call invert( mat2 )

                if ( k == kto ) then
                  ! Linear interpolation.
                  !call lin_sc( dtm( f ), psim, psi0 )
                  IQUV( :, f ) = matmul( mat2, matmul( mat1, IQUV( :, f ) ) + psim * sm( :, f ) + psi0 * s0( :, f ) )
                else
                  ! Parabolic interpolation.
                  call par_sc( dtm( f ), dtp( f ), psim, psi0, psip )
                  IQUV( :, f ) = matmul( mat2, matmul( mat1, IQUV( :, f ) ) + psim * sm( :, f ) + psi0 * s0( :, f ) + psip * sp( :, f ) )
                endif

              enddo ! 1 <= f <= spectrum_size
            end if ! is_horizontal_ray

            ! Increment J00 and J20
            continue

          enddo ! kfrom <= k <= fto, with kstep, z-index of layers

          ! Find maximal optical depth in each multiplet.
          do line = 1, atom%ntran
            begin_f = multiplets( line )%begin
            end_f   = multiplets( line )%end
            tau( line ) = maxval( tau_nu( begin_f:end_f ) )
          end do ! 1 <= line <= atom%ntran
          write(*, '(a, *(f5.2, x))') 'max tau(:) = ', tau

          ! Normalize IQUV
          !if ( mu > horizontal_mu_limit ) then
          !  write(*, *) angle, mu, maxval( IQUV( 1, 1:200 ) )
          !  IQUV( 2, 1:200 ) = IQUV( 2, 1:200 ) / maxval( IQUV( 1, 1:200 ) )
          !  IQUV( 3, 1:200 ) = IQUV( 3, 1:200 ) / maxval( IQUV( 1, 1:200 ) )
          !  IQUV( 4, 1:200 ) = IQUV( 4, 1:200 ) / maxval( IQUV( 1, 1:200 ) )
          !  IQUV( 1, 1:200 ) = IQUV( 1, 1:200 ) / maxval( IQUV( 1, 1:200 ) )
          !  IQUV( 1:4, 1:200 ) = IQUV( 1:4, 1:200 ) * 1d2
            call gp%options( 'set style data linespoints; set grid; set key top left; set colorsequence classic' )
            call gp%options( 'set xrange[-1.1:2.2];' )
            write(tmp_str, '(i4)') angle
            call gp%title( 'angle = ' // tmp_str  )
            call gp%plot( multiplets(1)%d_lambdas, tau_nu( 1:200 ) )
          !  call gp%plot( multiplets(1)%d_lambdas, transpose( IQUV( 1:4, 1:200 ) ) )
          !  !call gp%title( 'I/I_max'  )
          !  !call gp%title( 'Q/I_max'  )
          !  !call gp%plot( multiplets(1)%d_lambdas, IQUV( 2, 1:200 ) )
          !  !call gp%title( 'U/I_max'  )
          !  !call gp%plot( multiplets(1)%d_lambdas, IQUV( 3, 1:200 ) )
          !  !call gp%title( 'V/I_max'  )
          !  !call gp%plot( multiplets(1)%d_lambdas, IQUV( 4, 1:200 ) )
          !endif

        enddo ! angle (RT FS)
        !=======================================================================

      stop
      enddo ! iteration & tolerance


        do while (iteration < 50 .and. maxval( tolerance ) > 1d-3 )

            J00 = 0.d0
            J20 = 0.d0
            J00_nu = 0.d0
            J20_nu = 0.d0

! Solve the RT equation and calculate the tensors J00 and J20
            do i = 1, in_fixed%no
                call calculate_tensors(slab, i, J00_nu(:,i), J20_nu(:,i), in_fixed%gammad)
            enddo

! Carry out the integration over frequency weighted by the line profile                
            do i = 1, slab%n_layers
                J00(i) = J00(i) + int_tabulated(in_observation%freq, J00_nu(i,:)*prof)
                J20(i) = J20(i) + int_tabulated(in_observation%freq, J20_nu(i,:)*prof)
            enddo

! Put the new values of nbar and omega
            slab%nbar(:,1) = J00 * (in_fixed%wl*1.d-8)**3 / (2.d0*PC*PH)
            slab%omega(:,1) = J20 / J00 * sqrt(2.d0)

            iteration = iteration + 1

            !+DEBUG
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', iteration - 1, ', nbar:     ', slab%nbar
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', iteration - 1, ', nbar_old: ', slab%nbar_old
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', iteration - 1, ', omega:     ', slab%omega
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', iteration - 1, ', omega_old: ', slab%omega_old
            write (*, '(a, i4, a, es20.8, a, i4)') 'Iteration: ', iteration - 1, ', max |dn/n|: ', maxval(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar)), ' at layer ', maxloc(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar), 1)
            write (*, '(a, i4, a, es20.8, a, i4)') 'Iteration: ', iteration - 1, ', max |do/o|: ', maxval(abs(slab%omega - slab%omega_old) / abs(slab%omega)), ' at layer ', maxloc(abs(slab%omega - slab%omega_old) / abs(slab%omega), 1)
            !-DEBUG

            tolerance(1) = maxval(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar))
            tolerance(2) = maxval(abs(slab%omega - slab%omega_old) / abs(slab%omega))

            slab%nbar_old = slab%nbar
            slab%omega_old = slab%omega

            !print *, 'Maximum relative change : ', tolerance
            !write (*, '(a, i4, a, 2es9.2)') 'Iteration: ', iteration - 1, ', max. rel. change: ', tolerance

        enddo

! Recompute the propagation matrix and emission vector for the last time
! for the synthesis of the Stokes profiles

        do loop_shell = 1, slab%n_layers

            ! For the final formal solution, take the LoS projection of the
            ! vertical velocity in the current slab.  Hazel has a stupid
            ! historical convention (like many other things in observational
            ! astronomy) that vmacro is positive for redshifted profiles.  This
            ! is opposite to the positive direction of the angles and the
            ! magnetic field as the positive vmacro points downwards, against
            ! the OZ direction.  So, negate the projection to conform with this
            ! convention.
            in_params%vmacro = -slab%v_z( loop_shell ) * cos( in_fixed%thetad * PI / 180d0 )

            ! To fill and solve the statistical equilibrium equations, restore 
            ! the magnetic field vector in the current slab and...
            in_params%bgauss  = slab%B(    loop_shell )
            in_params%thetabd = slab%thB(  loop_shell )
            in_params%chibd   = slab%chiB( loop_shell )

            ! ...pass the obtained nbar and omega in the current shell.
            nbarExternal  = slab%nbar(  loop_shell, : )
            omegaExternal = slab%omega( loop_shell, : )

            ! For the first component only...
            call fill_SEE( in_params, in_fixed, 1, error )

! Calculate the absorption/emission coefficients for a given transition
! TODO : set velocity to zero for velocity-free approximation
            call calc_rt_coef(in_params, in_fixed, in_observation, 1)
                        
            if (.not.allocated(epsI)) allocate(epsI(in_fixed%no))
            if (.not.allocated(epsQ)) allocate(epsQ(in_fixed%no))
            if (.not.allocated(epsU)) allocate(epsU(in_fixed%no))
            if (.not.allocated(epsV)) allocate(epsV(in_fixed%no))
            if (.not.allocated(etaI)) allocate(etaI(in_fixed%no))
            if (.not.allocated(etaQ)) allocate(etaQ(in_fixed%no))
            if (.not.allocated(etaU)) allocate(etaU(in_fixed%no))
            if (.not.allocated(etaV)) allocate(etaV(in_fixed%no))
            if (.not.allocated(rhoQ)) allocate(rhoQ(in_fixed%no))
            if (.not.allocated(rhoU)) allocate(rhoU(in_fixed%no))
            if (.not.allocated(rhoV)) allocate(rhoV(in_fixed%no))
            if (.not.allocated(dtau)) allocate(dtau(in_fixed%no))
            

            if (in_fixed%use_atomic_pol == 1) then
! Emission                
                epsI = epsilon(0,:)
                epsQ = epsilon(1,:)
                epsU = epsilon(2,:)
                epsV = epsilon(3,:)
                
! Absorption including stimulated emission
                etaI = eta(0,:) - use_stim_emission_RT * eta_stim(0,:)
                etaQ = eta(1,:) - use_stim_emission_RT * eta_stim(1,:)
                etaU = eta(2,:) - use_stim_emission_RT * eta_stim(2,:)
                etaV = eta(3,:) - use_stim_emission_RT * eta_stim(3,:)

! Magneto-optical effects
                if (use_mag_opt_RT == 1) then
                    rhoQ = mag_opt(1,:) - use_stim_emission_RT * mag_opt_stim(1,:)
                    rhoU = mag_opt(2,:) - use_stim_emission_RT * mag_opt_stim(2,:)
                    rhoV = mag_opt(3,:) - use_stim_emission_RT * mag_opt_stim(3,:)
                else
                    rhoQ = 0.d0
                    rhoU = 0.d0
                    rhoV = 0.d0
                endif
            else
! Emission
                epsI = epsilon_zeeman(0,:)
                epsQ = epsilon_zeeman(1,:)
                epsU = epsilon_zeeman(2,:)
                epsV = epsilon_zeeman(3,:)

! Absorption including stimulated emission
                etaI = eta_zeeman(0,:) - use_stim_emission_RT * eta_stim_zeeman(0,:) + 1.d-20
                etaQ = eta_zeeman(1,:) - use_stim_emission_RT * eta_stim_zeeman(1,:) + 1.d-20
                etaU = eta_zeeman(2,:) - use_stim_emission_RT * eta_stim_zeeman(2,:) + 1.d-20
                etaV = eta_zeeman(3,:) - use_stim_emission_RT * eta_stim_zeeman(3,:) + 1.d-20

! Magneto-optical terms
                if (use_mag_opt_RT == 1) then
                    rhoQ = mag_opt_zeeman(1,:) - use_stim_emission_RT * mag_opt_stim_zeeman(1,:)
                    rhoU = mag_opt_zeeman(2,:) - use_stim_emission_RT * mag_opt_stim_zeeman(2,:)
                    rhoV = mag_opt_zeeman(3,:) - use_stim_emission_RT * mag_opt_stim_zeeman(3,:)
                else
                    rhoQ = 0.d0
                    rhoU = 0.d0
                    rhoV = 0.d0
                endif
            endif

! Set emission vector at this shell
            slab%emission_vector(1,loop_shell,:, angle) = epsI
            slab%emission_vector(2,loop_shell,:, angle) = epsQ
            slab%emission_vector(3,loop_shell,:, angle) = epsU
            slab%emission_vector(4,loop_shell,:, angle) = epsV

! Set propagation matrix at this shell
            do i = 1, in_fixed%no
                call fill_absorption_matrix(slab%propagation_matrix(:,:,loop_shell,i, angle),&
                    etaI(i),etaQ(i),etaU(i),etaV(i),rhoQ(i),rhoU(i),rhoV(i))
            enddo

! Multiply by the density
            slab%propagation_matrix(:,:,loop_shell,:, angle) = slab%propagation_matrix(:,:,loop_shell,:, angle) * &
                slab%density(loop_shell)
            slab%emission_vector(:,loop_shell,:, angle) = slab%emission_vector(:,loop_shell,:,angle) * &
                slab%density(loop_shell)

        enddo

        call synthesize_stokes(slab, in_params, in_fixed, output)

        open(unit=18,file='Jbar_tensors.dat',action='write',status='replace')
        write(18,*) slab%n_layers
        do i = 1, slab%n_layers
            !write(18,*) maxval(slab%tau(i,:)), slab%nbar(i,1), slab%omega(i,1)
            !+DEBUG
            !write (*, *) 'max tau, nbar, omega: ', maxval( slab%tau(i,:) ), slab%nbar(i,1), slab%omega(i,1)
            !-DEBUG
        enddo
        close(18)

        deallocate(J00_nu)
        deallocate(J20_nu)
        deallocate(J00)
        deallocate(J20)
        deallocate(prof)

      deallocate( slab%propagation_matrix, slab%emission_vector, IQUV )


    end subroutine do_transfer

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------    
    subroutine calculate_tensors(slab, freq, J00, J20, gamma)
    type(type_slab) :: slab
    integer :: freq
    real(kind=8) :: formal_sol_polarized(4), Inten(4), J00(:), J20(:), gamma
    integer :: k, km, kp, loop_mu, loop_direction
    real(kind=8) :: chim, chi0, chip, dtp, dtm, exu
    real(kind=8) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp, mu, Qtilde
    
    integer :: i, j, n, nmus, kfrom, kto, kstep, dummy
    real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:)
    real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4)
    
        n = slab%n_layers
        !nmus = slab%nmus
        nmus = slab%aq_size
        
        allocate(ab_matrix(4,4,n))
        allocate(source_vector(4,n))
                                
! Transform K into K* and then into K'
        do i = 1, 4
            do j = 1, 4
                ab_matrix(i,j,:) = slab%propagation_matrix( i, j, :, freq, dummy ) / slab%propagation_matrix( 1, 1, :, freq, dummy )
            enddo
            ab_matrix( i, i, : ) = ab_matrix( i, i, : ) - 1.d0
            source_vector(i,:) = slab%emission_vector( i, :, freq, dummy ) / slab%propagation_matrix(1,1,:,freq, dummy)
        enddo

        J00 = 0.d0
        J20 = 0.d0

        do loop_mu = 1, nmus

            !mu = slab%mus(loop_mu)
            mu = cos( slab%aq_inclination( loop_mu ) * PI / 180d0 )

! If going upwards
            if (mu > 0) then
                kfrom = 2
                kto = n
                kstep = 1
            else
! If going downwards
                kfrom = n-1
                kto = 1
                kstep = -1
            endif

! Boundary condition
            Inten = slab%boundary(:,loop_mu,freq)

! Calculate J00
            J00(kfrom-kstep) = J00(kfrom-kstep) + J00_FACTOR * slab%aq_weight(loop_mu) * Inten(1)

! Calculate J20 taking into account the contribution of I, Q and U
            Qtilde = cos(2.d0*gamma*PI/180.d0) * Inten(2) - sin(2.d0*gamma*PI/180.d0) * Inten(3)
            J20(kfrom-kstep) = J20(kfrom-kstep) + J20_FACTOR * slab%aq_weight(loop_mu) * &
                ( (3.d0*mu**2-1.d0) * Inten(1) - 3.d0*(1.d0-mu**2) * Qtilde )

            do k = kfrom, kto

! Parabolic short-characteristics
                if (k /= kto) then
                    km = k - 1
                    kp = k + 1
                    chim = slab%propagation_matrix(1,1,km,freq, dummy)                    
                    chi0 = slab%propagation_matrix(1,1,k,freq, dummy)                    
                    chip = slab%propagation_matrix(1,1,kp,freq, dummy)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq, dummy)
                    chi0 = slab%propagation_matrix(1,1,k,freq, dummy)
                    chip = 0.d0
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = 0.d0
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = 0.d0
                endif
                    
                dtm = 0.5d0 * (chim + chi0) * dm
                dtp = 0.5d0 * (chi0 + chip) * dp
                    
                if (dtm >= 1.d-4) then
                    exu = dexp(-dtm)
                else
                    exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
                endif
                
                call lin_sc(dtm,psim_lin,psi0_lin)
                mat1 = exu * identity_4x4 - psim_lin*ab_matrix(:,:,km)
                mat2 = identity_4x4 + psi0_lin * ab_matrix(:,:,k)
                call invert(mat2)
                        
                if (k /= kto) then
                    call par_sc(dtm,dtp,psim,psi0,psip)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
                else
                    call lin_sc(dtm,psim,psi0)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
                endif

! Calculate J00
                J00(k) = J00(k) + J00_FACTOR * slab%aq_weight(loop_mu) * Inten(1)

! Calculate J20 taking into account the contribution of I, Q and U
                Qtilde = cos(2.d0*gamma*PI/180.d0) * Inten(2) - sin(2.d0*gamma*PI/180.d0) * Inten(3)
                J20(k) = J20(k) + J20_FACTOR * slab%aq_weight(loop_mu) * &
                    ( (3.d0*mu**2-1.d0) * Inten(1) - 3.d0*(1.d0-mu**2) * Qtilde )
                        
            enddo
        enddo

        deallocate(ab_matrix)
        deallocate(source_vector)

    end subroutine calculate_tensors

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------    
    subroutine synthesize_stokes(slab, in_params, in_fixed, output)
    type(type_slab) :: slab
    type(variable_parameters) :: in_params
    type(fixed_parameters) :: in_fixed
    real(kind=8) :: output(0:3,in_fixed%no)
    real(kind=8) :: formal_sol_polarized(4), Inten(4)
    integer :: k, km, kp, freq
    real(kind=8) :: chim, chi0, chip, dtp, dtm, exu
    real(kind=8) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp, mu
    
    integer :: i, j, n, nmus, kfrom, kto, kstep, dummy
    real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:), total_tau(:)
    real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4), Ic
    
        n = slab%n_layers
        !nmus = slab%nmus
        
        allocate(ab_matrix(4,4,n))
        allocate(source_vector(4,n))
        allocate(total_tau(in_fixed%no))

        total_tau = 0.d0
        !tau = 0.d0

        do freq = 1, in_fixed%no
                                
! Transform K into K* and then into K'
            do i = 1, 4
                do j = 1, 4
                    ab_matrix(i,j,:) = slab%propagation_matrix(i,j,:,freq, dummy) / slab%propagation_matrix(1,1,:,freq, dummy)
                enddo
                ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
                source_vector(i,:) = slab%emission_vector(i,:,freq, dummy) / slab%propagation_matrix(1,1,:,freq, dummy)
            enddo

            ! The mu direction cosine points to the observer at (thetad, chid, gammad).
            mu = cos( in_fixed%thetad * PI / 180d0 )

! If going upwards
            kfrom = 2
            kto = n
            kstep = 1

            ! Stokes vector at the lower boundary: no polarization, intensity is
            ! from Eq. (12.33) from Landi degl'Innocenti & Landolfi (2004).
            ! The last formula is expanded for simplicity.
            Inten(1)   = I0_allen( in_fixed%wl, sqrt( 1d0 + ( mu**2 - 1d0 ) * ( 1d0 + in_params%height / RSUN )**2 ) )
            Inten(2:4) = 0d0

            do k = kfrom, kto

! Parabolic short-characteristics
                if (k /= kto) then
                    km = k - 1
                    kp = k + 1
                    chim = slab%propagation_matrix(1,1,km,freq, dummy)
                    chi0 = slab%propagation_matrix(1,1,k,freq, dummy)                    
                    chip = slab%propagation_matrix(1,1,kp,freq, dummy)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq, dummy)
                    chi0 = slab%propagation_matrix(1,1,k,freq, dummy)
                    chip = 0.d0
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = 0.d0
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = 0.d0
                endif
                    
                dtm = 0.5d0 * (chim + chi0) * dm
                dtp = 0.5d0 * (chi0 + chip) * dp

                total_tau(freq) = total_tau(freq) + dtm
                !slab%tau(k,freq) = slab%tau(km,freq) + dtm
                    
                if (dtm >= 1.d-4) then
                    exu = dexp(-dtm)
                else
                    exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
                endif
                
                call lin_sc(dtm,psim_lin,psi0_lin)
                mat1 = exu * identity_4x4 - psim_lin*ab_matrix(:,:,km)
                mat2 = identity_4x4 + psi0_lin * ab_matrix(:,:,k)
                call invert(mat2)
                        
                if (k /= kto) then
                    call par_sc(dtm,dtp,psim,psi0,psip)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
                else
                    call lin_sc(dtm,psim,psi0)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
                endif
                        
            enddo

            output(0:3,freq) = Inten
            
        enddo

        print *, 'Maximum tau=', maxval(total_tau)
        open(unit=18,file='final_tau.dat',action='write',status='replace')
        write(18,*) maxval(total_tau)
        write(18,*) total_tau
        close(18)

! Return I/Ic, Q/Ic, U/Ic and V/Ic
        Ic = output(0,1)
        do i = 0, 3
            output(i,:) = output(i,:) / Ic
        enddo

        deallocate(ab_matrix)
        deallocate(source_vector)

    end subroutine synthesize_stokes

end module nonlinear_transfer
