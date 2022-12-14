module all
    use iso_fortran_env
    implicit none

    private
    public :: everything
    type everything

        integer :: nx, ny, Tinfo = 540, Tfile = 349, niter = 0, piter = 0
        real(kind=real64) :: dx, dy, h, norm_res, err = 1.e-2, k = 1.0, Pr = 0.1
        real(kind=real64) :: Ra, a_dif, a_adv, total_time, dt_dif, dt_adv, dt, time = 0., size,  pi = 4.D0*DATAN(1.D0)
        real(kind=real64), allocatable :: T(:, :), S(:, :), W(:, :), vx(:, :), vy(:, :), x(:)
        real(kind=real64), dimension(1, 1) :: Vmax
        character(len = 50):: Tinit = 'none'
        character(len = 100) :: fname = 'lowPrandt.txt'

    contains

        procedure, public :: read_inputs, Initialize, VelocitySolve, Timestep, WriteTField, WriteInfos
        procedure, private :: setBC_T, setBC_0, ComputeVx, ComputeVy, &
            Vgrad, Laplacian, gradTx

        end type everything
        
contains

    ! ===== Initialize =====

        subroutine read_inputs(a)
            implicit none
            class(everything) :: a
            integer :: nx, ny
            real(kind=real64) :: Ra, Pr, a_dif, a_adv, err, total_time

            character(len=50) :: Tinit
            ! import parameter from text file
            namelist /inputs/ Pr, nx, ny, total_time, Ra, err, a_dif, a_adv, Tinit
            if (command_argument_count() > 0) &
                call get_command_argument(1, a%fname)

            ! open text files store the data
            open(1, file = a%fname, status = 'old')
            read(1, inputs)
            close(1)

            a%Pr = Pr
            a%nx = nx
            a%ny = ny
            a%Ra = Ra
            a%total_time = total_time
            a%a_dif = a_dif
            a%a_adv = a_adv
            a%err = err
            a%Tinit = Tinit

        end subroutine  

        subroutine initialize(a)
            implicit none

            class(everything) :: a
            integer :: i = 0, j = 0

            allocate(a%T(a%nx, a%ny), a%S(a%nx, a%ny), a%W(a%nx, a%ny), a%vx(a%nx, a%ny), a%vy(a%nx, a%ny), a%x(a%nx))

            a%T = 0; a%S = 0; a%W = 0; a%vx = 0; a%vy = 0

            ! Open file to store T
            if (a%Tinit == 'random') then
                open(a%Tfile, file = 'results/Tfield_random.txt')
                open(a%Tinfo, file = 'results/info_random.csv')
            else if (a%Tinit == 'cosine') then
                open(a%Tfile, file = 'results/Tfield_cosine.txt')
                open(a%Tinfo, file = 'results/info_cosine.csv')
            end if

            ! Initialize dimensions 
            a%h = 1./(a%ny - 1)
            a%dx = a%h
            a%dy = a%h
            a%dt_dif = a%a_dif * a%h**2 / a%k

            ! Create x vector to generate initial model
            a%x(1) = 0.
            do i = 2, a%nx
                a%x(i) = a%x(i-1) + a%h
            end do
            a%x(a%nx) = 1.

            ! Generate initial model
            a%size = a%nx * a%h
            if (a%Tinit == 'random') then
                call random_number(a%T)
            else if (a%Tinit == 'cosine') then
                do concurrent(i = 1:a%nx, j = 1:a%ny)
                    a%T(i, j) = 1./2 * (1 + cos(3*a%pi*a%x(i)/a%size))
                end do
            end if

            ! Set BC to 1 on the bottom, 0 on the top and on the boarders.
            call a%setBC_T()

            ! Initialize vorticity
            ! call a%gradTx(a%T, a%grad_Tx)
            a%W = -a%dt * a%Pr * a%Ra * gradTx(a, a%T)
            call a%SetBC_0(a%W)

            ! First T field update
            a%T = a%T + a%dt*Laplacian(a, a%T)
            call a%SetBC_T()    
            
            call a%WriteTField()

        end subroutine initialize

    ! ===== Linear algebra =====
    ! ----- Public ----- (mid level)

        subroutine VelocitySolve(a)

            implicit none
            class(everything) :: a
            
            ! poisson solver to get S from W
            a%norm_res = Vcycle_2DPoisson(a%S, a%W, a%h)
            do while (a%err <= a%norm_res)
                a%norm_res = Vcycle_2DPoisson(a%S, a%W, a%h)
            end do
            call a%SetBC_0(a%S)
            call a%ComputeVx(a%S, a%vx)
            call a%ComputeVy(a%S, a%vy)

            a%Vmax(1, 1) = maxval(sqrt(a%vx**2 + a%vy**2))

        end subroutine VelocitySolve

        subroutine Timestep(a)
            implicit none
            class(everything) :: a

            ! calculate advective timestep and min(adv, diff)
            a%dt_adv = a%a_adv*min(a%dx/maxval(a%vx), a%dy/maxval(a%vy))
            a%dt = min(a%dt_dif, a%dt_adv)

            ! update vorticity
            a%W = a%W + a%dt * (a%Pr * Laplacian(a, a%W) - Vgrad(a, a%W) &
                - a%Pr * a%Ra * gradTx(a, a%T))
            call a%SetBC_0(a%W)

            ! take advection and diffusion time step
            a%T = a%T + a%dt*(a%k * Laplacian(a, a%T) - Vgrad(a, a%T))
            call a%setBC_T()
        end subroutine Timestep
    
    ! ----- Private -----

        function Laplacian(a, input)

            implicit none
            class(everything) :: a
            integer :: i, j
            real(kind=real64) :: input(a%nx, a%ny), Laplacian(a%nx, a%ny)

            do i = 1, a%nx
                do j = 1, a%ny
                    ! Set BC to 0
                    if (i == 1 .or. i == a%nx .or. j == 1 .or. j == a%ny) then
                        Laplacian(i, j) = 0
                    else ! Compute second derivative
                        Laplacian(i, j) = (input(i+1, j) - 2*input(i, j) + input(i-1, j))/a%dx**2 + &
                            (input(i, j-1) - 2*input(i, j) + input(i, j+1))/a%dy**2
                    end if
                end do
            end do

        end function Laplacian

        function gradTx(a, input)

            implicit none
            class(everything) :: a
            integer :: i, j
            real(kind=real64) :: input(a%nx, a%ny), gradTx(a%nx, a%ny)

            do j = 1, a%ny
                do i = 1, a%nx
                    if (i == a%nx) then
                        gradTx(i, j) = 0
                    else
                        gradTx(i, j) = (input(i+1, j) - input(i-1, j))/(a%dx*2)
                    end if
                end do
            end do

        end function gradTx

        subroutine ComputeVx(a, input, output)

            implicit none
            class(everything) :: a
            integer :: i, j
            real(kind=real64) :: input(a%nx, a%ny), output(a%nx, a%ny)
            
            ! Calculate center derivative of a 2D array (in this cas vy on dy, using stream function Phi)
            do i = 1, a%nx
                do j = 1, a%ny
                    if (i == 1 .or. i == a%nx .or. j == 1 .or. j == a%ny) then
                        output(i, j) = 0
                    else
                        output(i, j) = (input(i, j+1) - input(i, j-1))/(2*a%dy)
                    end if
                end do
            end do

        end subroutine ComputeVx
        
        subroutine ComputeVy(a, input, output)

            implicit none
            class(everything) :: a
            integer :: i, j
            real(kind=real64) :: input(a%nx, a%ny), output(a%nx, a%ny)
            
            ! Calculate center derivative of a 2D array (in this cas vy on dy, using stream function Phi)
            do i = 1, a%nx
                do j = 1, a%ny
                    if (i == 1 .or. i == a%nx .or. j == 1 .or. j == a%ny) then
                        output(i, j) = 0
                    else
                        output(i, j) = (input(i-1, j) - input(i+1, j))/(2*a%dx)
                    end if
                end do
            end do

        end subroutine ComputeVy
        
        function Vgrad(a, input)

            implicit none
            class(everything) :: a
            integer :: i, j
            real(kind=real64) :: input(a%nx, a%ny), Vgrad(a%nx, a%ny), &
                UpVx(a%nx, a%ny), UpVy(a%nx, a%ny)
            
            ! Calculate v * grad T or W
            do i = 1, a%nx
                do j = 1, a%ny
                    if (i == 1 .or. i == a%nx .or. j == 1 .or. j == a%ny) then
                        UpVx(i, j) = 0
                        UpVy(i, j) = 0
                    else
                        ! Calculate upwind derivative for Vx
                        if (a%vx(i, j) > 0) then
                            UpVx(i, j) = a%vx(i, j) * (input(i, j) - input(i-1, j))/a%dx
                        else
                            UpVx(i, j) = a%vx(i, j) * (input(i+1, j) - input(i, j))/a%dx
                        end if

                        ! Calculate upwind derivative for Vy
                        if (a%vy(i, j) > 0) then
                            UpVy(i, j) = a%vy(i, j) * (input(i, j) - input(i, j-1))/a%dy
                        else
                            UpVy(i, j) = a%vy(i, j) * (input(i, j+1) - input(i, j))/a%dy
                        end if
                    end if
                end do
            end do
            Vgrad = UpVx + UpVy

        end function Vgrad
    
    ! ===== Velocity solver ======
    
        recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
            ! Here we add the given recursive function using V cycles
            implicit none
            real*8 resV
            real*8,intent(inout):: u_f(:,:)  ! arguments
            real*8,intent(in)   :: rhs(:,:),h
            integer         :: nx,ny,nxc,nyc, i  ! local variables
            real*8,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
            real*8            :: alpha=0.7, res_rms
            nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
            if( nx-1/=2*((nx-1)/2) .or. ny-1/=2*((ny-1)/2) ) &
            stop 'ERROR:not a power of 2'
            nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size
            
            if (min(nx,ny)>5) then  ! not the coarsest level
                allocate(res_f(nx,ny),corr_f(nx,ny), &
                corr_c(nxc,nyc),res_c(nxc,nyc))
                !---------- take 2 iterations on the fine grid--------------
                res_rms = iteration_2DPoisson(u_f,rhs,h,alpha) 
                res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
                !---------- restrict the residue to the coarse grid --------
                call residue_2DPoisson(u_f,rhs,h,res_f) 
                call restrict(res_f,res_c)
                !---------- solve for the coarse grid correction -----------
                corr_c = 0.  
                res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2) ! *RECURSIVE CALL
                !---- prolongate (interpolate) the correction to the fine grid 
                call prolongate(corr_c,corr_f)
                !---------- correct the fine-grid solution -----------------
                u_f = u_f - corr_f  
                !---------- two more smoothing iterations on the fine grid---
                res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
                res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
                deallocate(res_f,corr_f,res_c,corr_c)
            else  
                !----- coarsest level (ny=5): iterate to get 'exact' solution
                do i = 1,100
                    res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
                end do
            end if

            resV = res_rms   ! returns the rms. residue

        end function Vcycle_2DPoisson

        function iteration_2DPoisson(u, f, h, alpha)
            implicit none

            ! does one iteration on u field
            ! is a function that returns the rms residue 

            real*8, intent(inout) :: u(:, :)
            real*8, intent(in) :: f(:, :), h, alpha
            real*8 :: R, sum_r, sum_f, iteration_2DPoisson, rms_f, rms_r
            integer :: nx, ny, i, j
            nx = size(f, 1); ny = size(f, 2)

            sum_r = 0
            sum_f = 0
            do concurrent (i = 1:nx, j = 1:ny)
                if (i==1 .or. i==nx .or. j==1 .or. j==ny) then 
                    R = 0
                else
                    R = (u(i, j+1) + u(i, j-1) + u(i+1, j) + u(i-1, j) - 4*u(i, j))/h**2 - f(i, j)
                    u(i, j) = u(i, j) + alpha * R * h**2 / 4
                end if
                sum_f = f(i, j)**2 + sum_f
                sum_r = R**2 + sum_r
            end do

            rms_f = sqrt(sum_f/(nx*ny))
            rms_r = sqrt(sum_r/(nx*ny))

            iteration_2DPoisson = rms_r/rms_f
                
        end function
        
        subroutine residue_2DPoisson(u, f, h, res)
            implicit none

            ! calculate the residue in array res

            real*8, intent(out) ::res(:, :)
            real*8, intent(in) ::  u(:, :), f(:, :), h
            integer :: nx, ny, i, j

            nx = size(f, 1); ny = size(f, 2)
            
            do concurrent (i = 1:nx, j = 1:ny)
                if (i==1 .or. i==nx .or. j==1 .or. j==ny) then 
                    res(i, j) = 0
                else
                    res(i, j) = (u(i, j+1) + u(i, j-1) + u(i+1, j) + u(i-1, j) - 4*u(i, j))/h**2 - f(i, j)
                end if
            end do
                    
        end subroutine
        
        subroutine restrict(fine, coarse)        
            ! copies every other points in fine into coarse grid
            
            real*8, intent(in) :: fine(:, :)
            real*8, intent(out) :: coarse(:, :)
            
            integer :: nx, ny, i, j, is, js
            
            nx = size(fine, 1); ny = size(fine, 2)
            is = 1; js = 1
            
            do i = 1, nx, 2
                do j = 1, ny, 2
                    coarse(is, js) = fine(i, j)
                    js = js + 1
                end do
                js = 1
                is = is + 1
            end do
                    
        end subroutine
        
        subroutine prolongate(coarse, fine)        
            ! copies coarse into every other point in fine
            ! does linear interpolation to fill the other points
            real*8, intent(out) :: fine(:, :)
            real*8, intent(in) :: coarse(:, :)
            
            integer :: nx, ny, i, j, nxx, nyy

            nx = size(coarse, 1); ny = size(coarse, 2)
            nxx = 2*nx - 1
            nyy = 2*ny - 1

            do concurrent (i = 1:nxx:2, j = 1:nyy:2)
                fine(i, j) = coarse((i+1)/2, (j+1)/2)
            end do

            ! interpolation
            do concurrent (i = 2:nxx-1:2, j = 1:nyy:2)
                fine(i, j) = (fine(i+1, j) + fine(i-1, j))/2.
            end do

            do concurrent (i = 1:nxx, j = 2:nyy-1:2)
                fine(i, j) = (fine(i, j+1) + fine(i, j-1))/2.
            end do

        end subroutine

    ! ===== Utilities =====

        subroutine setBC_T(a)
        
            implicit none
            class(everything) :: a
        
            a%T(:, 1) = 1.
            a%T(:, a%ny) = 0
            a%T(1, :) = a%T(2, :)
            a%T(a%nx, :) = a%T(a%nx-1, :)
        
        end subroutine setBC_T
        
        subroutine setBC_0(a, input)
        
            implicit none
            class(everything) :: a
            real(kind=real64) :: input(a%nx, a%ny)
        
            input(:, 1) = 0
            input(:, a%ny) = 0
            input(1, :) = 0
            input(a%nx, :) = 0
        
        end subroutine setBC_0
        
        subroutine WriteTField(a)

            ! Write a 2D array into a text file
            implicit none
            class(everything) :: a
            integer :: j = 0

            do j = 1, a%ny ! Don't forget to change the maximum nx/ny number !
                write(a%Tfile,'(100000(1pe13.5))') a%T(:, j)
            end do

        end subroutine
        
        subroutine WriteInfos(a)

            implicit none
            class(everything) :: a

            write(a%Tinfo,'(*(G0.7,:,","))') 'Tinit', 'Pr', 'Ttime', 'nx', 'ny'
            write(a%Tinfo,'(*(G0.7,:,","))') a%Tinit, a%Pr, a%total_time, a%nx, a%ny
            
        end subroutine
    
end module all