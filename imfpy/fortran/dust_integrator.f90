module dust_integrator_mod
    use iso_fortran_env, only : dp => real64
    implicit none
contains

    subroutine integrate_particles(num_particles, num_steps, dt, gm_sun, q_over_m, beta, e_field, b_field, initial_state, results, status)
        integer, intent(in) :: num_particles
        integer, intent(in) :: num_steps
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: gm_sun
        real(dp), intent(in) :: q_over_m(num_particles)
        real(dp), intent(in) :: beta(num_particles)
        real(dp), intent(in) :: e_field(3)
        real(dp), intent(in) :: b_field(3)
        real(dp), intent(in) :: initial_state(6, num_particles)
        real(dp), intent(out) :: results(6, num_particles, num_steps)
        integer, intent(out) :: status

        integer :: p, step
        real(dp) :: state(6)
        real(dp) :: k1(6), k2(6), k3(6), k4(6)
        real(dp) :: temp_state(6)

        if (num_particles <= 0 .or. num_steps <= 1 .or. dt <= 0.0_dp) then
            status = -1
            results = 0.0_dp
            return
        end if

        status = 0

!$omp parallel do default(shared) private(p, step, state, k1, k2, k3, k4, temp_state)
        do p = 1, num_particles
            state(:) = initial_state(:, p)
            results(:, p, 1) = state(:)

            do step = 2, num_steps
                call compute_derivative(state, gm_sun, q_over_m(p), beta(p), e_field, b_field, k1)

                temp_state = state + 0.5_dp * dt * k1
                call compute_derivative(temp_state, gm_sun, q_over_m(p), beta(p), e_field, b_field, k2)

                temp_state = state + 0.5_dp * dt * k2
                call compute_derivative(temp_state, gm_sun, q_over_m(p), beta(p), e_field, b_field, k3)

                temp_state = state + dt * k3
                call compute_derivative(temp_state, gm_sun, q_over_m(p), beta(p), e_field, b_field, k4)

                state = state + (dt / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
                results(:, p, step) = state
            end do
        end do
!$omp end parallel do

    end subroutine integrate_particles

    subroutine compute_derivative(state, gm_sun, qom, beta, e_field, b_field, deriv)
        real(dp), intent(in) :: state(6)
        real(dp), intent(in) :: gm_sun
        real(dp), intent(in) :: qom
        real(dp), intent(in) :: beta
        real(dp), intent(in) :: e_field(3)
        real(dp), intent(in) :: b_field(3)
        real(dp), intent(out) :: deriv(6)

        real(dp) :: position(3)
        real(dp) :: velocity(3)
        real(dp) :: r, inv_r3, grav_coeff
        real(dp) :: acceleration(3)
        real(dp) :: lorentz(3)
        real(dp) :: cross(3)
        real(dp), parameter :: tiny = 1.0e-12_dp

        position = state(1:3)
        velocity = state(4:6)

        r = sqrt(position(1)**2 + position(2)**2 + position(3)**2)
        inv_r3 = 0.0_dp
        if (r > tiny) then
            inv_r3 = 1.0_dp / (r**3)
        end if

        grav_coeff = -gm_sun * (1.0_dp - beta) * inv_r3
        acceleration = grav_coeff * position

        cross(1) = velocity(2) * b_field(3) - velocity(3) * b_field(2)
        cross(2) = velocity(3) * b_field(1) - velocity(1) * b_field(3)
        cross(3) = velocity(1) * b_field(2) - velocity(2) * b_field(1)
        lorentz = qom * (e_field + cross)

        deriv(1:3) = velocity
        deriv(4:6) = acceleration + lorentz
    end subroutine compute_derivative

    subroutine version_string(buffer)
        character(len=*), intent(out) :: buffer
        character(len=*), parameter :: version = "dust_integrator/1.0"

        if (len(buffer) < len(version)) then
            buffer = version(1:len(buffer))
        else
            buffer(1:len(version)) = version
            if (len(buffer) > len(version)) then
                buffer(len(version)+1:) = ' '
            end if
        end if
    end subroutine version_string

end module dust_integrator_mod
