module nonlinear_transfer
use vars
use SEE
use rt_coef
use maths
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
    integer :: i, loop_iteration, loop_shell, error, j
    real(kind=8) :: I0, Q0, U0, V0, ds, Imax, Ic, factor, eta0, psim, psi0
    real(kind=8) :: dnum, relative_change(2), vmacro
    real(kind=8), allocatable, dimension(:) :: epsI, epsQ, epsU, epsV, etaI, etaQ, etaU, etaV, dtau
    real(kind=8), allocatable, dimension(:) :: rhoQ, rhoU, rhoV, delta, prof(:)
    real(kind=8), allocatable :: StokesM(:), kappa_prime(:,:), kappa_star(:,:), identity(:,:)
    real(kind=8), allocatable :: source(:), m1(:,:), m2(:,:), Stokes0(:)
    real(kind=8), allocatable :: O_evol(:,:), psi_matrix(:,:), J00(:), J20(:), J00_nu(:,:), J20_nu(:,:)

    integer(kind=4) :: line, layer, angle
    real(kind=8) :: mu, illumination_cone_cosine, illumination_cone_sine, cos_alpha

      if ( verbose_mode == 1 ) then
          print *, 'Starting transfer...'
      endif

      ! Create arrays for the current and the previous radiation field
      ! parameters, $\bar{n}$ and $\omega$, ...
      allocate( slab%nbar(      slab%n_layers, atom%ntran ) )
      allocate( slab%omega(     slab%n_layers, atom%ntran ) )
      allocate( slab%nbar_old(  slab%n_layers, atom%ntran ) )
      allocate( slab%omega_old( slab%n_layers, atom%ntran ) )

      ! ...and initialize the current nbar and omega in all layers with the same
      ! values using Allen's tables to start iterating.
      if ( verbose_mode == 1 ) then
          write (*, '(a)') 'Initial nbar and omega:'
      endif
      do line = 1, atom%ntran
        slab%nbar(  :, line ) = nbar_allen(  atom%wavelength( line ), in_fixed, in_params, atom%reduction_factor(       line ) )
        slab%omega( :, line ) = omega_allen( atom%wavelength( line ), in_fixed, in_params, atom%reduction_factor_omega( line ) )
        if ( verbose_mode == 1 ) then
            write (*, '(4x, a, f9.3, a, e8.2, a, e8.2)') 'wavelength = ', atom%wavelength( line ), ' A, nbar = ', slab%nbar( 1, line ), ', omega = ', slab%omega( 1, line )
        endif
      enddo
      ! Set the previous values to the current ones.
      slab%nbar_old  = slab%nbar
      slab%omega_old = slab%omega

      ! Create an array for the boundary conditions...
      allocate( slab%boundary( 4, slab%aq_size, in_fixed%no ), source = 0d0 )

      ! The illumination cone for a given height has
      illumination_cone_sine   = RSUN / ( RSUN + in_params%height )
      illumination_cone_cosine = sqrt( 1d0 - illumination_cone_sine**2 )
      ! from Eq. (12.32) of Landi degl'Innocenti & Landolfi (2004).

      ! ...and use Allen's tables for the center-to-limb variation of the
      ! continuum intensity to set all Stokes I values for the current angle.
      if ( verbose_mode == 1 ) then
          write (*, '(a)') 'Boundary intensities:'
      endif
      do angle = 1, slab%aq_size
        ! for all inclinations $\theta$ having $\mu = \cos\theta > \cos\gamma$
        mu = cos( slab%aq_inclination( angle ) )
        if ( mu > illumination_cone_cosine ) then
          ! Eq. 12.33 of Landi degl'Innocenti & Landolfi (2004).
          cos_alpha = sqrt( mu**2 - illumination_cone_cosine**2 ) / illumination_cone_sine
          slab%boundary( 1, angle, : ) = I0_allen( in_fixed%wl, cos_alpha )
        else
          ! For other angles there is no incoming boundary intensity and it is
          ! already set to zero.
          continue
        endif
        if ( verbose_mode == 1 ) then
          write (*, '(4x, a, f6.3, a, e8.2)') 'mu = ', mu, ', I0 = ', slab%boundary( 1, angle, 1 )
        endif
      enddo
      stop


      allocate(slab%emission_vector(4,slab%n_layers,in_fixed%no))
      allocate(slab%propagation_matrix(4,4,slab%n_layers,in_fixed%no))

      allocate(slab%tau(slab%n_layers,in_fixed%no))

      allocate(prof(in_fixed%no))

      allocate(J00_nu(slab%n_layers,in_fixed%no))
      allocate(J20_nu(slab%n_layers,in_fixed%no))
      allocate(J00(slab%n_layers))
      allocate(J20(slab%n_layers))

      relative_change = 100.d0

      loop_iteration = 1

! velocity-free approximation
        vmacro = in_params%vmacro
        in_params%vmacro = 0.0
        
        do while (loop_iteration < 50 .and. relative_change(1) > 1d-3 .and. relative_change(2) > 1d-3)

            do loop_shell = 1, slab%n_layers

                ! To fill and solve the statistical equilibrium equations,
                ! restore  the magnetic field vector in the current slab and...
                in_params%bgauss  = slab%B(    loop_shell )
                in_params%thetabd = slab%thB(  loop_shell )
                in_params%chibd   = slab%chiB( loop_shell )

                ! ...pass the obtained nbar and omega in the current shell.
                nbarExternal  = slab%nbar(  loop_shell, : )
                omegaExternal = slab%omega( loop_shell, : )

                ! For the first component only...
                call fill_SEE( in_params, in_fixed, 1, error )

                ! To calculate the absorption/emission coefficients, restore the
                ! damping parameters...
                in_params%damping = slab%damping(  loop_shell )
                in_params%vdopp   = slab%vthermal( loop_shell )
                ! TODO: The LoS velocity is needed as well but it should be made
                ! with a different angle quadrature.
                !in_params%vmacro  = slab%vmacro( loop_shell )

                call calc_rt_coef( in_params, in_fixed, in_observation, 1 )

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
            
                if (.not.allocated(StokesM)) allocate(StokesM(4))

                if (.not.allocated(source)) allocate(source(4))
                if (.not.allocated(kappa_star)) allocate(kappa_star(4,4))

                StokesM(1) = in_fixed%Stokes_incident(0)
                StokesM(2) = in_fixed%Stokes_incident(1)
                StokesM(3) = in_fixed%Stokes_incident(2)
                StokesM(4) = in_fixed%Stokes_incident(3)

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
                slab%emission_vector(1,loop_shell,:) = epsI
                slab%emission_vector(2,loop_shell,:) = epsQ
                slab%emission_vector(3,loop_shell,:) = epsU
                slab%emission_vector(4,loop_shell,:) = epsV

! Set propagation matrix at this shell
                do i = 1, in_fixed%no
                    call fill_absorption_matrix(slab%propagation_matrix(:,:,loop_shell,i),&
                        etaI(i),etaQ(i),etaU(i),etaV(i),rhoQ(i),rhoU(i),rhoV(i))
                enddo

! Multiply by the density
                slab%propagation_matrix(:,:,loop_shell,:) = slab%propagation_matrix(:,:,loop_shell,:) * &
                    slab%density(loop_shell)
                slab%emission_vector(:,loop_shell,:) = slab%emission_vector(:,loop_shell,:) * &
                    slab%density(loop_shell)

            enddo

            J00 = 0.d0
            J20 = 0.d0
            J00_nu = 0.d0
            J20_nu = 0.d0

! Generate line profile
            dnum = in_params%vdopp*1.d5 / (in_fixed%wl*1.d-8)
            prof = profile(in_params%damping,in_observation%freq / dnum) / (dnum*SQRTPI)

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

            loop_iteration = loop_iteration + 1

            !+DEBUG
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', loop_iteration - 1, ', nbar:     ', slab%nbar
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', loop_iteration - 1, ', nbar_old: ', slab%nbar_old
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', loop_iteration - 1, ', omega:     ', slab%omega
            !write (*, '(a, i4, a, 2es20.8)') 'Iteration: ', loop_iteration - 1, ', omega_old: ', slab%omega_old
            write (*, '(a, i4, a, es20.8, a, i4)') 'Iteration: ', loop_iteration - 1, ', max |dn/n|: ', maxval(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar)), ' at layer ', maxloc(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar), 1)
            write (*, '(a, i4, a, es20.8, a, i4)') 'Iteration: ', loop_iteration - 1, ', max |do/o|: ', maxval(abs(slab%omega - slab%omega_old) / abs(slab%omega)), ' at layer ', maxloc(abs(slab%omega - slab%omega_old) / abs(slab%omega), 1)
            !-DEBUG

            relative_change(1) = maxval(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar))
            relative_change(2) = maxval(abs(slab%omega - slab%omega_old) / abs(slab%omega))

            slab%nbar_old = slab%nbar
            slab%omega_old = slab%omega

            !print *, 'Maximum relative change : ', relative_change
            !write (*, '(a, i4, a, 2es9.2)') 'Iteration: ', loop_iteration - 1, ', max. rel. change: ', relative_change

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
            
            if (.not.allocated(StokesM)) allocate(StokesM(4))

            if (.not.allocated(source)) allocate(source(4))
            if (.not.allocated(kappa_star)) allocate(kappa_star(4,4))

            StokesM(1) = in_fixed%Stokes_incident(0)
            StokesM(2) = in_fixed%Stokes_incident(1)
            StokesM(3) = in_fixed%Stokes_incident(2)
            StokesM(4) = in_fixed%Stokes_incident(3)

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
            slab%emission_vector(1,loop_shell,:) = epsI
            slab%emission_vector(2,loop_shell,:) = epsQ
            slab%emission_vector(3,loop_shell,:) = epsU
            slab%emission_vector(4,loop_shell,:) = epsV

! Set propagation matrix at this shell
            do i = 1, in_fixed%no
                call fill_absorption_matrix(slab%propagation_matrix(:,:,loop_shell,i),&
                    etaI(i),etaQ(i),etaU(i),etaV(i),rhoQ(i),rhoU(i),rhoV(i))
            enddo

! Multiply by the density
            slab%propagation_matrix(:,:,loop_shell,:) = slab%propagation_matrix(:,:,loop_shell,:) * &
                slab%density(loop_shell)
            slab%emission_vector(:,loop_shell,:) = slab%emission_vector(:,loop_shell,:) * &
                slab%density(loop_shell)

        enddo

        call synthesize_stokes(slab, in_params, in_fixed, output)

        open(unit=18,file='Jbar_tensors.dat',action='write',status='replace')
        write(18,*) slab%n_layers
        do i = 1, slab%n_layers
            write(18,*) maxval(slab%tau(i,:)), slab%nbar(i,1), slab%omega(i,1)
            !+DEBUG
            write (*, *) 'max tau, nbar, omega: ', maxval( slab%tau(i,:) ), slab%nbar(i,1), slab%omega(i,1)
            !-DEBUG
        enddo
        close(18)

        deallocate(J00_nu)
        deallocate(J20_nu)
        deallocate(J00)
        deallocate(J20)
        deallocate(prof)
            
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
    
    integer :: i, j, n, nmus, kfrom, kto, kstep
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
                ab_matrix(i,j,:) = slab%propagation_matrix(i,j,:,freq) / slab%propagation_matrix(1,1,:,freq)
            enddo
            ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
            source_vector(i,:) = slab%emission_vector(i,:,freq) / slab%propagation_matrix(1,1,:,freq)
        enddo

        J00 = 0.d0
        J20 = 0.d0

        do loop_mu = 1, nmus

            !mu = slab%mus(loop_mu)
            mu = cos( slab%aq_inclination( loop_mu ) )

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
                    chim = slab%propagation_matrix(1,1,km,freq)                    
                    chi0 = slab%propagation_matrix(1,1,k,freq)                    
                    chip = slab%propagation_matrix(1,1,kp,freq)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)
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
    
    integer :: i, j, n, nmus, kfrom, kto, kstep
    real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:), total_tau(:)
    real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4), Ic
    
        n = slab%n_layers
        !nmus = slab%nmus
        
        allocate(ab_matrix(4,4,n))
        allocate(source_vector(4,n))
        allocate(total_tau(in_fixed%no))

        total_tau = 0.d0
        slab%tau = 0.d0

        do freq = 1, in_fixed%no
                                
! Transform K into K* and then into K'
            do i = 1, 4
                do j = 1, 4
                    ab_matrix(i,j,:) = slab%propagation_matrix(i,j,:,freq) / slab%propagation_matrix(1,1,:,freq)
                enddo
                ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
                source_vector(i,:) = slab%emission_vector(i,:,freq) / slab%propagation_matrix(1,1,:,freq)
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
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)                    
                    chip = slab%propagation_matrix(1,1,kp,freq)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)
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
                slab%tau(k,freq) = slab%tau(km,freq) + dtm
                    
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
